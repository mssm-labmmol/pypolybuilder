import unittest
import sys
sys.path.append('/home/vitor/dendribuilder')
import utils
from atom import Atom
from bond import Bond
from angle import Angle


class bondTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
    	bondTests.a1 = Atom(1,'H',None,None,None,None,None,None)
    	bondTests.a2 = Atom(2,'NT',None,None,None,None,None,None)
        bondTests.a3 = Atom(3,'NT',None,None,None,None,None,None)

    def testFindBondParam(self):
    	bond = Bond(bondTests.a1,bondTests.a2,None,None)
    	self.assertNotEqual(bond.get_param,None)

    def testContainsAtom(self):
        bond = Bond(bondTests.a1,bondTests.a2,None,None)
        self.assertTrue(bond.contains_atom(bondTests.a1))
        self.assertTrue(bond.contains_atom(bondTests.a2))

    def testNotContainsAtom(self):
        bond = Bond(bondTests.a1,bondTests.a2,None,None)
        self.assertFalse(bond.contains_atom(bondTests.a3))

    def testGetOtherAtom(self):
        bond = Bond(bondTests.a1,bondTests.a2,None,None)
        self.assertEqual(bond.get_other_atom(bondTests.a1),bondTests.a2)

    def testWrongGetOtherAtom(self):
        bond = Bond(bondTests.a1,bondTests.a2,None,None)
        self.assertNotEqual(bond.get_other_atom(bondTests.a1),bondTests.a1)