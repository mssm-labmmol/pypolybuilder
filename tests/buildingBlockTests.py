import unittest
import sys
sys.path.append('/home/vitor/dendribuilder')
import utils
from atom import Atom
from bond import Bond
from angle import Angle
from BuildingBlock import BuildingBlock

class buildingBlockTests(unittest.TestCase):

    def testFindBondParam(self):
    	bb = BuildingBlock('bb_oh.itp')
        self.assertNotEqual(bb.get_atom_list(),None)
        self.assertNotEqual(bb.get_bond_list(),None)