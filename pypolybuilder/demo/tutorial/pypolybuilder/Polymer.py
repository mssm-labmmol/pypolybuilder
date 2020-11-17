'''
Created on Dec 22, 2013

@author: root
'''

from Itp      import Itp
from Angle    import Angle
from Bond     import Bond
from Dihedral import Dihedral
import random
import copy
import Utils
from Utils import ngen
from Utils import connections_paths
from Utils import bbs_paths
from BuildingBlockFactory import BuildingBlockFactory
from Connection import Connection

class Polymer(Itp):
                                      
    def create_coords(self):
        for atom in self.get_atom_list():
            x = random.uniform(0, 1)
            y = random.uniform(0, 1)
            z = random.uniform(0, 1)
            atom.set_x(x)
            atom.set_y(y)
            atom.set_z(z)

    #Read connect file and stores connections and building blocks
    def read_connections(self):
    	f = open(Utils.connections_paths, 'r')
        line = " "
        for line in f:
        	if "[ BUILDING BLOCKS ]" in line:
        		for line in f:
        			if not line.isspace() and line[0] != ";"  and "[ CONNECTS ]" not in line:    
        				line = line.split()
        				bb_number  = int(line[0])
        				bb_resname = line[1]
        				building_block = self.bb_factories[bb_resname].makeBBForItp(self)
        				self.building_blocks[bb_number] = building_block
			    		self.connect_bb(building_block)
        				self.dict_number_resname[bb_number] = bb_resname
        			else:
        				break
        	if "[ CONNECTS ]" in line:
        		for line in f:
        			line = line.split()
        			connection = Connection(line[0],line[1],line[2],line[3])
        			self.connections.append(connection)
            
        f.close()

    #read and store the building blocks files paths
    def read_bbs(self):
    	for path in Utils.bbs_paths:
    		bb_factory = BuildingBlockFactory(path)
    		self.bb_factories[bb_factory.resName] = bb_factory

                
    def connect_bb(self,bb):
    	atoms     = copy.copy(bb.get_atom_list())
        bonds     = copy.copy(bb.get_bond_list())
        angles    = copy.copy(bb.get_angle_list())
        dihedrals = copy.copy(bb.get_dihedral_list())
        
        self.atom_list.extend(atoms)
        self.bond_list.extend(bonds)
        self.angle_list.extend(angles)
        self.dihedral_list.extend(dihedrals)

    def connect_bbs(self):
    	for connection in self.connections:    		    		
    		firstBB = self.building_blocks[connection.b1]
    		secondBB = self.building_blocks[connection.b2]

    		atom1 = firstBB.atom_list[connection.atom_b1-1]
    		atom2 = secondBB.atom_list[connection.atom_b2-1]

    		bond = Bond(atom1,atom2,2,None)
    		self.bond_list.append(bond)




    def __init__(self): 
    	self.connections  = []
    	self.bb_factories = {}
    	self.dict_number_resname = {}
    	self.atom_list      = []
    	self.bond_list      = []
    	self.angle_list     = []
    	self.dihedral_list  = []
    	self.exclusion_list = []
    	self.building_blocks = {}