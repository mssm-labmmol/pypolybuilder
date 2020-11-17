import unittest
import sys
sys.path.append('/home/vitor/dendribuilder')
import itp
import getopt
import copy
from polymer import Polymer
from dendrimer import Dendrimer
import utils

class PolymerTests(unittest.TestCase):
    def testReadBBs(self):
    	polymer = Polymer()
    	polymer.read_bbs()
    	self.assertNotEqual(polymer.bb_factories['H'],None)
    	self.assertNotEqual(polymer.bb_factories['CCO'],None)
    	self.assertNotEqual(polymer.bb_factories['OH'],None)

    def testReadConnections(self):
    	polymer = Polymer()
    	polymer.read_bbs()
    	polymer.read_connections()
    	self.assertEqual(len(polymer.connections),4)
