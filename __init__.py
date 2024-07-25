__version__ = 'v1.1'
__author__ = "Feng Group"
__email__ = ""

import warnings
warnings.simplefilter("ignore")
from ._cluster import cluster
from ._index import index
from ._strain import strain
from ._bootstrap import bootstrap
from ._species import species

