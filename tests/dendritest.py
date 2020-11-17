import unittest
import sys
sys.path.append('/home/vitor/dendribuilder')
import testEnv
import itp
import getopt
import copy
from polymer import Polymer
from dendrimer import Dendrimer
import utils
import polymerTests



def load_tests(loader, tests, pattern):
    print loader.discover('.')
    return loader.discover('.')

if __name__ == '__main__':
    all_tests = unittest.TestLoader().discover('.', pattern='*.py')
    unittest.TextTestRunner().run(all_tests)