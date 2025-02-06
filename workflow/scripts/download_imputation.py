import json
import requests
import datetime
import time
import logging
from urllib.parse import urlparse
import urllib
import sys
import os
import re

if 'snakemake' not in globals():
    import yaml
    with open("config/config.yaml", 'r') as ymlfile:
        cfg = yaml.safe_load(ymlfile)
    class snakemake_class_testing:
        input = []
        params = {}
    snakemake = snakemake_class_testing()
    snakemake.params['token'] = cfg["impute"]["imputation"]["default"]["token"]
    snakemake.params['outpath'] = "temp/sandbox"
    os.makedirs(snakemake.params['outpath'], exist_ok=True)
    snakemake.input = ["intermediate/imputation/imputation/nacc-gsa_EAS_imputation.json"]


# Configure logging
logging.basicConfig(
    format='%(asctime)s [%(levelname)s] %(message)s',
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)]
)


def format_speed(speed):
    """Makes download speed in human-readable."""
    units = ["B/s", "KB/s", "MB/s", "GB/s"]
    precision = 0
    while speed >= 1024 and len(units) > 1:
        speed /= 1024
        units.pop(0)
        precision += 1
    return f"{speed:.{min(precision, 3)}f} {units[0]}"
    

def progress(file_name, count, block_size, total_size):
    """Prints download progress every 15 seconds with file name."""
    global last_progress_time
    global start_time
    global reported

    now = time.time()

    if count == 0:
        last_progress_time = now
        start_time = now
        reported = False
        return
    size = int(count * block_size)
    progress_size = f"{size // (1024 * 1024)} MB"
    percent = min(int(count * block_size * 100 / total_size), 100)
    speed = format_speed(size / (now - start_time))

    initial_report = False
    if now - start_time >= 2 and not reported:
        initial_report = True
        reported = True

    if now - last_progress_time >= 15 or initial_report:  # Update every 15 seconds
        last_progress_time = now
        logging.info(f"Downloading {file_name}: {percent}% | {progress_size} | {speed}")

def jobinfo(settings, token=None):
    """Fetch job information from the API."""
    settings['token'] = token if token is not None else settings['token']
    r_check = requests.get(
        '{u}/jobs/{j}'.format(u=settings['url'], j=settings['id']),
        headers={'X-Auth-Token': settings['token']})
    if r_check.status_code != 200:
        token_config = snakemake.params['token']
        if r_check.status_code == 401:
            if token is not None or settings['token'] == token_config:
                raise ValueError('Bad or expired API token')
            else:
                return jobinfo(settings, token=token_config)
        else:
            raise Exception(f"GET /jobs/{settings['id']}/status {r_check.status_code}")
    jinfo = r_check.json()
    jinfo['settings'] = settings
    return jinfo

def fmt_delta(start_time):
    """Formats elapsed time as HH:MM:SS."""
    delt = datetime.datetime.now() - start_time
    hours = delt.seconds // 3600
    mins = (delt.seconds % 3600) // 60
    sec = (delt.seconds % 3600) % 60
    ts = f"{hours:02d}:{mins:02d}:{sec:02d}"
    return f"{delt.days} Days {ts}" if delt.days > 0 else ts

# Read job details
jsonfile = snakemake.input[0]
with open(jsonfile, 'r') as j:
    imputation = json.load(j)

jinfo = jobinfo(imputation)

if jinfo['state'] < 4:
    logging.info(f"Waiting for {imputation['id']} to complete.")

submission = datetime.datetime.fromtimestamp(jinfo['submittedOn'] / 1000)

friendly_jobname = re.sub(r'_submitted20\d\d-\d\d-\d\d\.\d+$', '', jinfo["name"])


# Monitor job status
while jinfo['state'] < 4:
    if jinfo['state'] == 1:
        print('{} waiting for {}. (Queue position {})'.format(
            imputation['id'], fmt_delta(submission), jinfo['positionInQueue']))
        time.sleep(60) # wait for 1 minute
    elif jinfo['state'] == 2:
        start_time = datetime.datetime.fromtimestamp(jinfo['startTime'] / 1000)
        logging.info(f"{imputation['id']} running for {fmt_delta(start_time)}.")
        time.sleep(20)
    elif jinfo['state'] == 3:
        logging.info(f"{imputation['id']} exporting.")
        time.sleep(20)
    else:
        raise Exception(f"{imputation['id']} died.")
    jinfo = jobinfo(jinfo['settings'])

# Prepare downloads
url_dl = 'https://{fqdn}/share/results/{{hash}}/{{name}}'.format(
    fqdn=urlparse(imputation['url']).netloc)

outputs = [{**x, 'download': [{'url_dl': url_dl.format(**y),
                               'filename': y['name'],
                               'size': y['size']} for y in x['files']]}
           for x in jinfo['outputParams']]

outputs_dict = {x['description']: x['download']
                for x in outputs if x['download']}

# download files
if jinfo['state'] not in [7, 10]:
    for desc, dl in outputs_dict.items():
        logging.info(f"Downloading {desc} for {friendly_jobname}")
        for file in dl:
            destfile = os.path.join(snakemake.params['outpath'], file['filename'])
            logging.info(f"Starting download: {file['filename']} for {friendly_jobname}")
            urllib.request.urlretrieve(
                file['url_dl'], filename=destfile,
                reporthook=lambda c, bs, ts: progress(f"{friendly_jobname} {file['filename']}", c, bs, ts))
        logging.info(f"Finished downloading {desc} for {friendly_jobname}")

# Save job info
with open(os.path.join(snakemake.params['outpath'], 'job.json'), 'w') as f:
    json.dump(jinfo, f)

# Handle job completion status
if jinfo['state'] in [4, 8]:
    logging.info(f"{imputation['id']} was successful.")
elif jinfo['state'] in [5, 6, 7, 9, 10]:
    state = {6: 'cancelled', 7: 'retired', 10: 'deleted',
             5: 'failed', 9: 'failed'}[jinfo['state']]
    state += '. Check logs' if jinfo['state'] in [5, 9] else '. Please rerun'
    logging.error(f"{imputation['id']} was {state}.")
    sys.exit(1)
else:
    logging.error(f"{imputation['id']} status not recognized. Status = {jinfo['state']}")
    sys.exit(1)
