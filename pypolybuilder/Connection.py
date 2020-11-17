
class Connection(object):
    def __init__(self, b1, b2, atom_b1, atom_b2, has_specified_dihedral=False, specified_dihedral=None):
        self.b1      = int(b1)
        self.b2      = int(b2)
        self.atom_b1 = int(atom_b1)
        self.atom_b2 = int(atom_b2)
        self.has_specified_dihedral = has_specified_dihedral
        self.specified_dihedral = specified_dihedral or {"atom_nrs": [], "funct": [], "ph0s": []}