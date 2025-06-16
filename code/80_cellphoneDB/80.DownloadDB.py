import pandas as pd
import glob
import os

from IPython.display import HTML, display
from cellphonedb.utils import db_releases_utils

cpdb_version = 'v4.1.0'

path_db = 'CellPhoneDB/db'
cpdb_target_dir = os.path.join(path_db, cpdb_version)

from cellphonedb.utils import db_utils
db_utils.download_database(cpdb_target_dir, cpdb_version)

