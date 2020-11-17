import unittest
import sys
sys.path.append('/home/vitor/dendribuilder')
import utils
from atom import Atom
from bond import Bond
from angle import Angle
from itp import Itp
from polymer import Polymer
from BuildingBlockFactory import BuildingBlockFactory

class buildingBlockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
    	buildingBlockTests.itp = Polymer()

    def testMakeBBForItp(self):
    	bbFactory = BuildingBlockFactory('bb_oh.itp')
    	bb = bbFactory.makeBBForItp(buildingBlockTests.itp)
    	self.assertNotEqual(bb,None)
    	self.assertEqual(bb.get_atom_list()[0].get_nr(),1)
    	buildingBlockTests.itp.set_atom_list(bb.get_atom_list())
    	bb2 = bbFactory.makeBBForItp(buildingBlockTests.itp)
        self.assertNotEqual(bb2,None)
    	self.assertEqual(bb2.get_atom_list()[0].get_nr(),3)
