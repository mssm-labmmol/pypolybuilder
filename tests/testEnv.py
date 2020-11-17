import unittest
import sys
sys.path.append('/home/vitor/dendribuilder')
from polymer import Polymer
from dendrimer import Dendrimer
import utils
from ff_params import ff_params

utils.connections_paths = 'connect.in'
utils.bbs_paths = ['bb_oh.itp','bb_cco.itp','bb_h.itp']
utils.ff_param = ff_params('list_param.itp')
