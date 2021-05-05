# This script will load the config from yaml
# Has to be in seperate modeul to avoid circular error

import sys
sys.path.append("python-packages/pyyaml")
import yaml

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

## Open the user master yaml file
with open('master_user_config.yaml') as f:

    data = yaml.load(f,Loader=Loader)
