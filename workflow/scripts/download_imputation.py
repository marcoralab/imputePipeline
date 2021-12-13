#!/usr/bin/env python

import json
import requests
import datetime
import time
from urllib.parse import urlparse
import urllib
import sys
import os

# jsonfile = 'intermediate/pre-impute_prep/imputation/EUR_imputation.json'
# imputation = {'url': 'https://imputationserver.sph.umich.edu/api/v2',
#               'token': 'eyJjdHkiOiJ0ZXh0XC9wbGFpbiIsImFsZyI6IkhTMjU2In0.eyJtYWlsIjoiZnVsdG9uaDFAZ21haWwuY29tIiwiZXhwaXJlIjoxNjQxMDUyNDc0NjMxLCJuYW1lIjoiQnJpYW4gRSBGdWx0b24tSG93YXJkIiwiYXBpIjp0cnVlLCJ1c2VybmFtZSI6ImJlZmgifQ.IhGvMBoHG12qCpjuAndPldZwGnRwBklB35xpQl2UjAM',
#               'id': 'job-20211207-143135-505'}
# imputation['id'] = 'job-20171123-122635-549'

def progress(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = min(int(count * block_size * 100 / total_size), 100)
    sys.stdout.write('\r  {}%, {} MB, {} KB/s, {} seconds passed'.format(
        percent, progress_size // (1024 * 1024), int(speed), int(duration)))


def jobinfo(settings, token=None):
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
            raise Exception('GET /jobs/{}/status {}'.format(
                settings['id'], r_check.status_code))
    jinfo = r_check.json()
    jinfo['settings'] = settings
    return jinfo


def fmt_delta(start_time):
    delt = datetime.datetime.now() - start_time
    hours = delt.seconds // 3600
    mins = (delt.seconds % 3600) // 60
    sec = (delt.seconds % 3600) % 60
    ts = '{h:02d}:{m:02d}:{s:02d}'.format(h=hours, m=mins, s=sec)
    if delt.days > 0:
        return '{D} Days {t}'.format(D=delt.days, t=ts)
    return ts


jsonfile = snakemake.input[0]
with open(jsonfile, 'r') as j:
    imputation = json.load(j)

jinfo = jobinfo(imputation)

if jinfo['state'] < 4:
    print('Waiting for {} to complete.\n\n'.format(imputation['id']))

submission = datetime.datetime.fromtimestamp(jinfo['submittedOn']/1000)

while jinfo['state'] < 4:
    print("\033[A                                      \033[A")
    if jinfo['state'] == 1:
        print('{} waiting for {}. (Queue position {})'.format(
            imputation['id'], fmt_delta(submission), jinfo['positionInQueue']))
        time.sleep(60) # wait for 1 minute
    elif jinfo['state'] == 2:
        start_time = datetime.datetime.fromtimestamp(jinfo['startTime']/1000)
        print('{} running for {}.'.format(
            imputation['id'], fmt_delta(start_time)))
        time.sleep(20) # wait for 20 sec
    elif jinfo['state'] == 2:
        print('{} exporting.'.format(imputation['id']))
        time.sleep(20) # wait for 20 sec
    else:
        raise Exception('{} died.'.format(imputation['id']))
    jinfo = jobinfo(jinfo['settings'])

url_dl = 'https://{fqdn}/share/results/{{hash}}/{{name}}'.format(
    fqdn=urlparse(imputation['url']).netloc)

outputs = [{**x, 'download': [{'url_dl': url_dl.format(**y),
                               'filename': y['name'],
                               'size': y['size']} for y in x['files']]}
           for x in jinfo['outputParams']]

outputs_dict = {x['description']: x['download']
                for x in outputs if x['download']}

if jinfo['state'] not in [7, 10]:
    for desc, dl in outputs_dict.items():
        print('Downloading ' + desc)
        for file in dl:
            destfile = os.path.join(snakemake.params['outpath'],
                                    file['filename'])
            print('\nFile: {} ...'.format(file['filename']))
            urllib.request.urlretrieve(
                file['url_dl'], filename=destfile, reporthook=progress)
        print('\n')

with open(os.path.join(snakemake.params['outpath'], 'job.json'), 'w') as f:
    json.dump(jinfo, f)

if jinfo['state'] in [4, 8]:
    print('{} was successful.'.format(imputation['id']))
elif jinfo['state'] in [5, 6, 7, 9, 10]:
    state = {6: 'cancelled', 7: 'retired', 10: 'deleted',
             5: 'failed', 9: 'failed'}[jinfo['state']]
    state += '. Check logs' if jinfo['state'] in [5, 9] else '. Please rerun'
    raise Exception('{} was {}.'.format(imputation['id'], state))
else:
    raise Exception('{} status not recognized. Status = {}'.format(
        imputation['id'], jinfo['state']))

def human_readable_size(size, decimal_places=4):
    for unit in ['bytes', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if size < 1024.0 or unit == 'PiB':
            break
        size /= 1024.0
    return f"{size:.{decimal_places}f} {unit}"

#human_readable_size(98558729)
