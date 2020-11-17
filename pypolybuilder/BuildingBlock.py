'''
Created on Apr 2, 2015

@author: vitor
'''
from Itp import Itp
from Atom     import Atom
from Bond     import Bond
from Angle    import Angle
from Dihedral import Dihedral
from Branch   import Branch
from Utils import DEBUG


class BuildingBlock(Itp):
    """
    classdocs
    """
    # Initializes the BuildingBlock setting his template
    def __init__(self, path):
        self.atom_list      = []
        self.bond_list      = []
        self.angle_list     = []
        self.dihedral_list  = []
        self.branch_list    = []
        self.exclusion_list = []
        self.resName        = ''
    
        with open(path, 'r') as f:     
            for line in f:
                if "[ moleculetype ]" in line:
                    self.read_resName(f)
                if "[ atoms ]" in line:
                    self.read_atoms(f)
                if "[ bonds ]" in line:
                    self.read_bonds(f)
                if "[ angles ]" in line:
                    self.read_angles(f)
                if "[ dihedrals ]" in line:
                    if DEBUG:
                        print("Found specifications for dihedrals in input blocks.")
                    self.read_dihedrals(f)
                if "[ branches ]" in line:
                    self.read_branches(f)
                if "[ exclusions ]" in line:
                    self.read_exclusions(f)

            for atom in self.get_atom_list():
                atom.set_bond_list([])
                atom.set_atom_list([])

            for bond in self.get_bond_list():
                bond.get_a_1().get_atom_list().append(bond.get_other_atom(bond.get_a_1()))
                bond.get_a_2().get_atom_list().append(bond.get_other_atom(bond.get_a_2()))
            # not sure, why we call this here:
            self.set_angle_list(self.create_all_angles(self.get_bond_list()))   # not sure, why we call this here
            self.set_dihedral_list(self.create_all_dihedrals())  # not sure, why we call this here

    def read_resName(self, f):
        for line in f:
            if line.isspace(): break
            if line[0] == ';': continue
            line = line.split()
            self.resName = line[0]

    def read_atoms(self, f):
        for line in f:
            if line.isspace(): break
            if line[0] == ';': continue
            line = line.split()
            atom = Atom(line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7])
            if len(line) > 8:
                if line[8] == "AR":
                    atom.set_aromatic()
            self.atom_list.append(atom)

    def read_bonds(self, f):
        for line in f:
            if line.isspace(): break
            if line[0] == ';': continue
            data = line.split()
            atom_1 = self.find_atom(data[0])
            atom_2 = self.find_atom(data[1])
            bond = Bond(atom_1,atom_2,data[2],data[3])
            self.bond_list.append(bond)

    def read_angles(self, f):
        for line in f:
            if line.isspace(): break
            if line[0] == ';': continue
            data = line.split()
            atom_1 = self.find_atom(data[0])
            atom_2 = self.find_atom(data[1])
            atom_3 = self.find_atom(data[2])
            angle = Angle(atom_1, atom_2, atom_3, data[3], data[4], True)
            self.angle_list.append(angle)
            
    def read_dihedrals(self, f):
        found_non_comments = False
        for line in f:
            if line.isspace():
                break
            if is_comment(line):
                continue
            data = line.split()
            atom_1 = self.find_atom(data[0])
            atom_2 = self.find_atom(data[1])
            atom_3 = self.find_atom(data[2])
            atom_4 = self.find_atom(data[3])
            dihedral = Dihedral(atom_1, atom_2, atom_3, atom_4, int(data[4]), data[5],
                                isDetermined=False, was_input=False)
            self.dihedral_list.append(dihedral)
            found_non_comments = True
        if found_non_comments:
            print("Dihedrals have been found in input. These will be printed in the Topology, but ignored for the "
                  "optimization of the torsions for the initial structure, eventually.")

    def read_exclusions(self, f):
        for line in f:
            if line.isspace(): break
            if line[0] == ';': continue
            data = line.split()
            atom_1 = self.find_atom(data[0])
            atom_2 = self.find_atom(data[1])
            atom_1.get_exclusion_list().append(int(atom_2.get_nr()))
            atom_1.get_exclusion_extra().append(int(atom_2.get_nr()))

    def read_branches(self,f):
        for line in f:
            if line.isspace(): break
            if line[0] == ';': continue
            line = line.split()
            branch = Branch(line[0],line[1])
            self.branch_list.append(branch)
                    
    def get_donor_branch(self):
        for branch in self.get_branch_list():
            if branch.get_donor() != 0:
                return branch
        
        return None


def is_comment(line, character=";"):
    return line.startswith(character)
