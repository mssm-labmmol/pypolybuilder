import unittest
import sys
sys.path.append('/home/vitor/dendribuilder')
import utils
from atom import Atom
from bond import Bond
from angle import Angle


class angleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
    	angleTests.a1 = Atom(1,'H',None,None,None,None,None,None)
    	angleTests.a2 = Atom(2,'NT',None,None,None,None,None,None)
    	angleTests.a3 = Atom(3,'H',None,None,None,None,None,None)
    	angleTests.a4 = Atom(4,None,None,None,None,None,None,None)
    	
    def testMakeFromBondsWithSuccess(self):    	
    	b1 = bond = Bond(angleTests.a1,angleTests.a2,2,None)
    	b2 = bond = Bond(angleTests.a2,angleTests.a3,2,None)
    	angle = Angle.makeFromBonds(b1,b2)
    	self.assertEqual(angle.get_a_1(),angleTests.a1)
    	self.assertEqual(angle.get_a_2(),angleTests.a2)
    	self.assertEqual(angle.get_a_3(),angleTests.a3)

    def testFindAngleParam(self):
    	angle = Angle(angleTests.a1,angleTests.a2,angleTests.a3,None,None)
    	self.assertNotEqual(angle.get_param,None)