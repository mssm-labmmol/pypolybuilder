'''
Created on Dec 22, 2013

@author: root
'''
import copy
import math
# from collections.Set import set
class Atom(object):
    __nr             = None
    __resnr          = None
    __cgnr           = None        
    __charge         = None
    __mass           = None
    __atomtype       = None
    __resid          = None
    __atom           = None
    __aromatic       = None
    __x              = None
    __y              = None
    __z              = None
    __bond_list      = None
    __atom_list      = None
    __exclusion_list = None
    __exclusion_extra= None
    __pairs_14_list = []
    __neighbor_list  = None
    __connected_list = None
    __is_in_specified_dihedral = None
    __specified_dihedral_vals = None
    __is_branchpoint = None
    __is_donor = None
    __is_acceptor = None
    __is_neighbor_of_tetraedral = None
    __was_central_atom = None
    __was_central_atom_for = None

    def __init__(self, nr, atomtype,resnr, resid,atom, cgnr, charge, mass, aromatic=False,
                 is_in_specified_dihedral=False, specified_dihedral_vals=None,
                 is_branchpoint=False, is_donor=False, is_acceptor=False, is_neighbor_of_tetraedral=False,
                 was_central_atom=False, was_central_atom_for=[]):
        self.__nr = nr
        self.__resnr = resnr
        self.__cgnr = cgnr
        self.__charge = charge
        self.__mass = mass
        self.__atomtype = atomtype
        self.__resid = resid
        self.__atom = atom
        self.__aromatic = aromatic
        self.__x    = None
        self.__y    = None
        self.__z    = None
        self.__bond_list = []
        self.__bond_list = copy.copy(self.__bond_list)
        self.__atom_list = []
        self.__atom_list = copy.copy(self.__atom_list)
        self.__exclusion_list = []
        self.__exclusion_extra = []
        self.__neighbor_list = []
        self.connected_list = set()
        self.__is_in_specified_dihedral = is_in_specified_dihedral
        self.__specified_dihedral_vals = specified_dihedral_vals or {"idx_in_dihed": None, "dihed_atom_nrs": [],
                                                                     "dihed_atoms": [], "funct": None, "ph0": "",
                                                                     "use_for_lookup": False}
        self.__is_branchpoint = is_branchpoint
        self.__is_donor = is_donor
        self.__is_acceptor = is_acceptor
        self.__is_neighbor_of_tetraedral = is_neighbor_of_tetraedral
        self.__was_central_atom = was_central_atom
        self.__was_central_atom_for = was_central_atom_for

    def get_pairs14_list(self):
        return self.__pairs_14_list

    def set_pairs14_list(self, value):
        pairs_14_list = copy.copy(value)
        for atom_nr in pairs_14_list:
            if (atom_nr <= int(self.get_nr())):
                value.remove(atom_nr)
        
        self.__pairs_14_list = value
    
    def get_exclusion_list(self):
        return self.__exclusion_list

    def get_exclusion_extra(self):
        return self.__exclusion_extra

    def set_exclusion_list(self, value):
        exclusion_list = copy.copy(value.get_exclusion_list())
        
        for atom_nr in exclusion_list:
            if(atom_nr <= int(self.get_nr())):
                value.get_exclusion_list().remove(atom_nr)
        
        self.__exclusion_list = value
        
    def set_exclusion_extra(self, value):
        exclusion_extra = copy.copy(value.get_exclusion_extra())
        
        for atom_nr in exclusion_extra:
            if(atom_nr <= int(self.get_nr())):
                value.get_exclusion_extra().remove(atom_nr)
        
        self.__exclusion_extra = value

    def del_exclusion_list(self):
        del self.__exclusion_list

    def get_atom_list(self):
        return self.__atom_list

    def get_atomtype(self):
        return self.__atomtype

    def get_resid(self):
        return self.__resid

    def get_is_in_specified_dihedral(self):
        return self.__is_in_specified_dihedral

    def get_specified_dihedral_vals(self):
        return self.__specified_dihedral_vals

    def is_aromatic(self):
        return self.__aromatic

    def set_aromatic(self):
        self.__aromatic = True

    def set_atom_list(self, value):
        self.__atom_list = value

    def set_was_central_atom(self, value):
        self.__was_central_atom = value

    def append_was_central_atom_for(self, value):
        self.__was_central_atom_for.append(value)

    def is_terminal(self):
        return len(self.__atom_list) < 2

    def get_is_branchpoint(self):
        return self.__is_branchpoint

    def get_is_acceptor(self):
        return self.__is_acceptor

    def get_is_donor(self):
        return self.__is_donor

    def del_atom_list(self):
        del self.__atom_list

    def get_was_central_atom(self):
        return self.__was_central_atom

    def get_was_central_atom_for(self):
        return self.__was_central_atom_for

    def get_bond_list(self):
        return self.__bond_list

    def set_bond_list(self, value):
        self.__bond_list = value

    def set_is_in_specified_dihedral(self, value):
        self.__is_in_specified_dihedral = value

    def set_specified_dihedral_vals(self, value):
        self.__specified_dihedral_vals = value

    def del_bond_list(self):
        del self.__bond_list

    def get_neighbor_list(self):
        return self.__neighbor_list

    def del_neighbor_list(self):
        del self.__neighbor_list

    def set_neighbor_list(self, value):
        self.__neighbor_list = value

    def set_is_branchpoint(self, value):
        self.__is_branchpoint = value

    def set_is_acceptor(self, value):
        self.set_is_branchpoint(True)
        self.__is_acceptor = value

    def set_is_donor(self, value):
        self.set_is_branchpoint(True)
        self.__is_donor = value

    def print_coordinates(self):
        print(str(self.get_x()) + " " + str(self.get_y()) + " " +str(self.get_z()))

    def print_bonds(self):
       print(len(self.get_bond_list()))
       #for bond in self.bond_list():
            #print(str(bond.get_a_1().get_nr()) +  " " + bond.get_a_1().get_atomtype() + " " + str(bond.get_a_2().get_nr()) + " " + bond.get_a_2().get_atomtype())

    def is_donor(self,branch_list):
        for branch in branch_list:
            if(branch.get_donor() == self.get_nr()):
                return True
        
        return False

        #Calculate the distance absolute value between two atoms
    def calcDist(self, atom2):
        xdiff = self.get_x() - atom2.get_x()
        ydiff = self.get_y() - atom2.get_y()
        zdiff = self.get_z() - atom2.get_z()

        distance = math.sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)

        return distance
    
    # Here we don't calculate the square root in order to save computer time 
    def CalcDistSquared(self, atom2):
        xdiff = self.get_x() - atom2.get_x()
        ydiff = self.get_y() - atom2.get_y()
        zdiff = self.get_z() - atom2.get_z()

        distsqr = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff

        return distsqr

    def calcNeighbor(self, atom2,cutoff2):
        if(self.CalcDistSquared(atom2) < cutoff2):
            self.addNeighbor(atom2)
            atom2.addNeighbor(self)
            
    def testAtomIsClose(self, atom2,cutoff2):
        if(self.CalcDistSquared(atom2) < cutoff2):
            return True

    def addNeighbor(self, atom):
        if(atom not in self.get_neighbor_list()):
            self.get_neighbor_list().append(atom)

    def setposition(self, x, y, z):
        # if (self.get_x() != None): return
        self.set_x(x)
        self.set_y(y)
        self.set_z(z)

    def get_nr(self):
        return int(self.__nr)

    def get_resnr(self):
        return self.__resnr

    def get_cgnr(self):
        return self.__cgnr
    
    def get_x(self):
        return self.__x

    def get_y(self):
        return self.__y

    def get_z(self):
        return self.__z

    def get_charge(self):
        return self.__charge

    def get_mass(self):
        return self.__mass

    def get_atomtype(self):
        return self.__atomtype

    def get_resid(self):
        return self.__resid

    def get_atom(self):
        return self.__atom

    def get_is_neighbor_of_tetraedral(self):
        return self.__is_neighbor_of_tetraedral

    def set_nr(self, value):
        self.__nr = value

    def set_resnr(self, value):
        self.__resnr = value

    def set_cgnr(self, value):
        self.__cgnr = value

    def set_charge(self, value):
        self.__charge = value

    def set_mass(self, value):
        self.__mass = value

    def set_is_neighbor_of_tetraedral(self, value):
        self.__is_neighbor_of_tetraedral = value

    def set_atomtype(self, value):
        self.__atomtype = value
        
    def set_x(self, value):
        self.__x = value

    def add_to_x(self,value):
        self.__x += value

    def set_y(self, value):
        self.__y = value

    def add_to_y(self, value):
        self.__y += value

    def set_z(self, value):
        self.__z = value

    def add_to_z(self, value):
        self.__z += value

    def set_resid(self, value):
        self.__resid = value


    def set_atom(self, value):
        self.__atom = value

    def __str__(self):
        return self.get_atomtype()

    def __repr__(self):
        return self.get_nr()

    nr = property(get_nr, set_nr, None, None)
    resnr = property(get_resnr, set_resnr, None, None)
    cgnr = property(get_cgnr, set_cgnr, None, None)
    charge = property(get_charge, set_charge, None, None)
    mass = property(get_mass, set_mass, None, None)
    atomtype = property(get_atomtype, set_atomtype, None, None)
    resid = property(get_resid, set_resid, None, None)
    atom = property(get_atom, set_atom, None, None)
    bond_list = property(get_bond_list, set_bond_list, del_bond_list, "bond_list's docstring")
    atom_list = property(get_atom_list, set_atom_list, del_atom_list, "atom_list's docstring")
    exclusion_list = property(get_exclusion_list, set_exclusion_list, del_exclusion_list, "exclusion_list's docstring")
    neighbor_list = property(get_neighbor_list, set_neighbor_list, del_neighbor_list, "neighbor_list's docstring")
    # connected_list = property(get_connected_list, set_connected_list, del_connected_list, "connected_list's docstring")


