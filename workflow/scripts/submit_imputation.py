import requests
import json
import time
import secrets
import string
import datetime
import logging
import sys

# Configure logging
logging.basicConfig(
    format='%(asctime)s [%(levelname)s] %(message)s',
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)]
)

if 'snakemake' not in globals():
    import yaml
    cohort = "colombian-psen1-e280a-mega_AMR"
    with open("config/config.yaml", 'r') as ymlfile:
        cfg = yaml.safe_load(ymlfile)["impute"]
    class snakemake_class_testing:
        input = {}
        params = {}
        wildcards = {}
        output = []
    imputation_defaults = {
        'nih': {
            'server': 'NIH',
            'refpanel': 'topmed-r3',
            'population': 'all'},
        'michigan': {
            'server': 'Michigan',
            'refpanel': 'hrc-r1.1',
            'population': 'mixed'}}
    if 'imputation' in cfg and 'default' in cfg['imputation']:
        default_cfg = cfg['imputation']['default']
        if ( 'server' in default_cfg and default_cfg['server'].lower() == 'michigan'
             or ('refpanel' in default_cfg
                 and default_cfg['refpanel'].lower() == 'hrc-r1.1')):
            default_imp = imputation_defaults['michigan']
        elif ('refpanel' in default_cfg and
              default_cfg['refpanel'].lower() == 'topmed-r3'):
            default_imp = imputation_defaults['nih']
        elif ('refpanel' in default_cfg and default_cfg['refpanel']
              and 'server' in default_cfg and default_cfg['server']
              and 'population' in default_cfg and default_cfg['population']):
            default_imp = default_cfg
        else:
            raise ValueError('Must specify at least server or panel')
        default_imp.update(default_cfg)
    else:
        default_imp = imputation_defaults['nih']
    imp_settings = default_imp.copy()
    if 'imputation' in cfg and cohort in cfg['imputation']:
        imp_settings.update(cfg['imputation'][cohort])
    if 'token' in imp_settings:
        token = imp_settings.pop('token')
    else:
        raise ValueError("Must provide either cohort or default API token.")
    snakemake = snakemake_class_testing()
    snakemake.params['server'] = imp_settings.pop('server')
    snakemake.params['token'] = token
    snakemake.params['imp_settings'] = imp_settings
    snakemake.input['contig_build'] = f"intermediate/imputation/rename_chrom/{cohort}_mapping.json"
    snakemake.input['vcf'] = [f"intermediate/imputation/ready/{cohort}_chr{x}.vcf.gz" for x in range(1, 23)]
    snakemake.wildcards['cohort'] = cohort
    snakemake.output.append(f"intermediate/imputation/{cohort}_imputation_new.json")

cohort = snakemake.wildcards['cohort']

def getjobs(url, token):
    r_jobs = requests.get(url + "/jobs", headers={'X-Auth-Token' : token })
    if r_jobs.status_code != 200:
        raise Exception('GET /jobs/ {}'.format(r_jobs.status_code))
    return r_jobs.json()['data']

if snakemake.params['server'].lower() == 'michigan':
    url = 'https://imputationserver.sph.umich.edu/api/v2'
elif snakemake.params['server'].lower() == 'nih':
    url = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
else:
    url = snakemake.params['server']

token = snakemake.params['token']

r = requests.get(url + "/jobs", headers={'X-Auth-Token': token })
if r.status_code == 401:
    raise ValueError('Bad or expired API token')

if r.status_code == 404:
    raise ValueError('Invalid Imputation Server API URL')
elif r.status_code != 200:
    raise Exception('Server Error: Status {}'.format(r_jobs.status_code))
else:
    try:
        r.json()
    except ValueError:
        raise ValueError('Invalid Imputation Server API URL or invalid response')

# get all jobs
jobs = getjobs(url, token)
incomplete = list(filter(lambda x: x['state'] < 4, jobs))

if len(incomplete) > 2:
    logging.info(f"{cohort} submission: 3 jobs already queued/running on server.")
    while len(incomplete) > 2:
        if max([x['state'] for x in incomplete]) == 1: #queued only
            qpos = min([x['positionInQueue'] for x in incomplete])
            logging.info(f"{cohort} submission: lowest queue position is {qpos}.")
        else:
            running = len([x for x in incomplete if x['state'] > 1])
            logging.info(f"{cohort} submission: {running} jobs running.")
        time.sleep(600) # wait for 10 minutes
        jobs = getjobs(url, token)
        incomplete = list(filter(lambda x: x['state'] < 4, jobs))
    logging.info(f"{cohort} submission: Job completed; ready to submit.")

# define password and job name, then remove extranious imputation params

data = snakemake.params['imp_settings'].copy()
data['password'] = ''.join(
    (secrets.choice(string.ascii_letters + string.digits)
     for i in range(48)))
data['job-name'] = '{}_submitted{}'.format(
    cohort,
    datetime.datetime.now().strftime("%Y-%m-%d.%H%M"))

if 'token' in data:
    del data['token']

if 'server' in data:
    del data['server']

# rudimentary build check

with open(snakemake.input['contig_build'], 'r') as f:
    contig_build = 'hg19' if json.load(f)['build'] == 37 else 'hg38'
    if 'build' in data and data['build'] != contig_build:
        raise ValueError('Specified build does not match contig names.')
    elif 'build' not in data and contig_build != 'hg19':
        raise ValueError('Build not specified but contig names not hg19.')

# submit new job

if url == 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2':
    submit = "/jobs/submit/imputationserver2"
else:
    submit = "/jobs/submit/minimac4"

r_submission = requests.post(url + submit,
    files=[('files', open(x, 'rb')) for x in snakemake.input['vcf']],
    data=data,
    headers={'X-Auth-Token': token })

json_submission = r_submission.json()
logging.info(f"{cohort} submission message: {json_submission['message']}")

if r_submission.status_code != 200:
    raise Exception('POST {} {}'.format(submit, r_submission.status_code))

# print id
logging.info(f"{cohort} submission id: {json_submission['id']}")

json_submission['password'] = data['password']
json_submission['job-name'] = data['job-name']
json_submission['url'] = url
json_submission['token'] = token
del data['password']
del data['job-name']
json_submission['settings'] = data

r_check = requests.get(
    '{u}/jobs/{j}/status'.format(u=url, j=json_submission['id']),
    headers={'X-Auth-Token': token })
if r_check.status_code != 200:
    raise Exception('GET /jobs/{}/status {}'.format(
        json_submission['id'], r_check.status_code))

logging.info(f"{cohort} queue position: {r_check.json()['positionInQueue']}")

serverstats = requests.get(url + '/server/counters').json()

logging.info(f"{cohort} queue status: "
             f"{serverstats['queue']['size']} jobs currently in queue.")

with open(snakemake.output[0], 'w') as jobinfo:
    json.dump(json_submission, jobinfo)
