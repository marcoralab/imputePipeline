import requests
import json
import time
import secrets
import string
import datetime

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
    print("Three jobs are already queued or running on the server:\n\n")
    while len(incomplete) > 2:
        print("\033[A                                      \033[A") #erase line
        if max([x['state'] for x in incomplete]) == 1: #queued only
            qpos = min([x['positionInQueue'] for x in incomplete])
            print("lowest queue position is {}.".format(qpos))
        else:
            running = len([x for x in incomplete if x['state'] > 1])
            print('{} jobs running.'.format(running))
        time.sleep(600) # wait for 10 minutes
        jobs = getjobs(url, token)
        incomplete = list(filter(lambda x: x['state'] < 4, jobs))
    print("Job completed; ready to submit.")

# define password and job name, then remove extranious imputation params

data=snakemake.params['imp_settings']
data['password'] = ''.join(
    (secrets.choice(string.ascii_letters + string.digits)
     for i in range(48)))
data['job-name'] = '{}_submitted{}'.format(
    snakemake.wildcards['cohort'],
    datetime.datetime.now().strftime("%Y-%m-%d.%H%M"))
if 'token' in data:
    del data['token']
if 'server' in data:
    del data['server']

# rudimentary build check

with open(snakemake.input.contig_build, 'r') as f:
    contig_build = 'hg19' if json.load(f)['build'] == 37 else 'hg38'
    if 'build' in data and data['build'] != contig_build:
        raise ValueError('Specified build does not match contig names.')
    elif 'build' not in data and contig_build != 'hg19':
        raise ValueError('Build not specified but contig names not hg19.')

# submit new job

if url == 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2':
    submit = "/jobs/submit/imputationserver"
else:
    submit = "/jobs/submit/minimac4"

r_submission = requests.post(url + submit,
    files=[('files', open(x, 'rb')) for x in snakemake.input.vcf],
    data=data,
    headers={'X-Auth-Token': token })
if r_submission.status_code != 200:
    print(r.json()['message'])
    raise Exception('POST {} {}'.format(submit, r_submission.status_code))

json_submission = r_submission.json()

# print message
print(json_submission['message'])
print(json_submission['id'])

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

print('Queue position : {}'.format(r_check.json()['positionInQueue']))

serverstats = requests.get(url + '/server/counters').json()

print('{} jobs currently running.'.format(serverstats['running']['runs']))

with open(snakemake.output[0], 'w') as jobinfo:
    json.dump(json_submission, jobinfo)
