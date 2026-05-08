import requests
import yaml

with open("config/config.yaml", 'r') as ymlfile:
    cfg = yaml.safe_load(ymlfile)["impute"]


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

def token_settings(token):
    params = {'cohorts': [k for k, v in cfg["imputation"].items() if v['token'] == token]}
    cohort = params['cohorts'][0]

    imp_settings = default_imp.copy()
    if cohort in cfg['imputation']:
        imp_settings.update(cfg['imputation'][cohort])
    if 'token' in imp_settings:
        token = imp_settings.pop('token')
    else:
        raise ValueError("Must provide either cohort or default API token.")

    params['server'] = imp_settings.pop('server')
    params['token'] = token
    params['imp_settings'] = imp_settings
    if params['server'].lower() == 'michigan':
        params['url'] = 'https://imputationserver.sph.umich.edu/api/v2'
    elif params['server'].lower() == 'nih':
        params['url'] = 'https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2'
    else:
        params['url'] = snakemake.params['server']
    return(params)

tokens =  {x['token'] for x in cfg["imputation"].values()}
accounts = {x: token_settings(x) for x in tokens}

i = 0
for k, v in accounts.items():
    accounts[k]["i"] = i
    i += 1
    accounts[k]['r'] = requests.get(v['url'] + "/jobs", headers={'X-Auth-Token': v['token'] })
    accounts[k]['status_code'] = accounts[k]['r'].status_code
    accounts[k]['response'] = None
    accounts[k]['fail'] = True
    if accounts[k]['status_code'] == 401:
        accounts[k]['problem'] = 'Bad or expired API token'
    elif accounts[k]['status_code'] == 404:
        accounts[k]['problem'] = 'Invalid Imputation Server API URL'
    elif accounts[k]['status_code'] != 200:
        accounts[k]['problem'] = 'Server Error: Status {}'.format(accounts[k]['status_code'])
    else:
        try:
            accounts[k]['response'] = accounts[k]['r'].json()
            accounts[k]['fail'] = False
            accounts[k]['problem'] = None
        except ValueError:
            accounts[k]['problem'] = 'Invalid Imputation Server API URL or invalid response'

output = [{'status_code': v['status_code'],
           'problem': v['problem'],
           'response': v['response'],
           'token': v['token'],
           'server': v['server'],
           'cohorts': v['cohorts']} for v in accounts.values()]

with open('server_status.yaml', 'w') as f:
    yaml.dump(output, f)