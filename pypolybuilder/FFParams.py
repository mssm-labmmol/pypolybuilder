"""
Created on Jan 22, 2015

@author: Vitor
"""
from BondParams import bond_param
from AngleParams import angle_param
from DihedralParams import dihe_param

COMMENT_CHARS = ["#"]


def quit_wrong_params(line_number, line):
    print("Wrong number of elements in line {}: \n{}".format(line_number, line))
    quit()


def strip_comments(line, comment_chars=COMMENT_CHARS):
        for comment_char in comment_chars:
            line = line.split(comment_char)[0]
        return line


def give_unique_identifier_to_dihedral_params(ff_param):
    # impropers
    dihed_params_list = ff_param.get_improper_params()
    for didx, dihed_param in enumerate(dihed_params_list):
        dihed_param.set_unique_identifier("imp{}".format(didx))
    ff_param.set_improper_params(dihed_params_list)
    # torsional dihedrals
    dihed_params_list = ff_param.get_dihe_params()
    for didx, dihed_param in enumerate(dihed_params_list):
        dihed_param.set_unique_identifier("tor{}".format(didx))
    ff_param.set_dihe_params(dihed_params_list)

    # elif dihe_key == "imp_branch":
    #     ff_params_list = Utils.ff_param.get_improper_branch_params()
    #     check_func = check_imp_dihed_atypes
    # elif dihe_key == "prop_branch":
    #     ff_params_list = Utils.ff_param.get_branch_diheparams()
    #     check_func = check_prop_dihed_atypes
    return ff_param


class ff_params(object):
    __bondparams  = None
    __angleparams = None
    __impparams   = None
    __diheparams  = None
    __transdihe   = None
    __cisdihe     = None
    __branch_diheparams = None
    __improper_branch_params = None
    __branch_diheds_given = None
    __branch_imps_given = None
    __double_params_given = False  # will be updated, if user has given fancy dihedrals parameters

    def __init__(self, params_filename):
        self.__bondparams  = []
        self.__angleparams = []
        self.__impparams = []
        self.__diheparams = []
        self.__trans = []
        self.__cis = []
        self.__branch_diheparams = []
        self.__improper_branch_params = []
        self.__branch_diheds_given = False
        self.__branch_imps_given = False
        if params_filename != '':
            self.read_bond_params(params_filename)
            self.read_conformations(params_filename)

    def get_bond_params(self):
        return self.__bondparams
    
    def get_angle_params(self):
        return self.__angleparams
        
    def get_improper_params(self):
        return self.__impparams
    
    def get_dihe_params(self):
        return self.__diheparams

    def get_branch_diheparams(self):
        return self.__branch_diheparams

    def get_improper_branch_params(self):
        return self.__improper_branch_params

    def get_branch_diheds_given(self):
        return self.__branch_diheds_given

    def get_branch_imps_given(self):
        return self.__branch_imps_given

    def get_double_params_given(self):
        return self.__double_params_given

    def set_double_params_given(self, value):
        self.__double_params_given = value

    def get_is_trans(self, dihedral):
        a1 = str(dihedral.get_a_1())
        a2 = str(dihedral.get_a_2())
        a3 = str(dihedral.get_a_3())
        a4 = str(dihedral.get_a_4())

        for ttype in self.__trans:
            if a1 == ttype[0] and a2 == ttype[1] and a3 == ttype[2] and a4 == ttype[3]:
                return True
        return False

    def get_is_cis(self, dihedral):
        a1 = str(dihedral.get_a_1())
        a2 = str(dihedral.get_a_2())
        a3 = str(dihedral.get_a_3())
        a4 = str(dihedral.get_a_4())

        for ttype in self.__cis:
            if a1 == ttype[0] and a2 == ttype[1] and a3 == ttype[2] and a4 == ttype[3]:
                return True
        return False

    def set_branch_diheds_given(self, value):
        print("Found specified dihedrals for dendrimer branch points!")
        self.__branch_diheds_given = value

    def set_branch_imps_given(self, value):
        print("Found specified improper dihedrals for dendrimer branch points!")
        self.__branch_imps_given = value

    def set_improper_params(self, value):
        self.__impparams = value

    def set_dihe_params(self, value):
        self.__diheparams = value

    def set_branch_diheparams(self, value):
        self.__branch_diheparams = value

    def set_improper_branch_params(self, value):
        self.__improper_branch_params = value

    def createConformation(self, line):
        line = line.split()
        a1 = line[0]
        a2 = line[1]
        a3 = line[2]
        a4 = line[3]
        conform = line[4]
        if conform == 'trans':
            self.__trans.append([a1, a2, a3, a4])
            return
        self.__cis.append([a1, a2, a3, a4])

    def read_conformations(self, filename):
        with open(filename, 'r') as f:
            reading_conformation = False
            for line in f:
                if reading_conformation:
                    if line.isspace():
                        break
                    self.createConformation(line)
                if 'CONFORMATION' in line:
                    reading_conformation = True

    def create_param(self, line_number, line, param_type):
        """THis functions creates parameters from columns in line"""
        elems = line.split()
        param = None
        if param_type == "bond":
            if len(elems) != 3:
                quit_wrong_params(line_number, line)
            return bond_param(*elems)  # atom_name_1, atom_name_2, param_name)
        elif param_type == "angle":
            if len(elems) != 4:
                quit_wrong_params(line_number, line)
            return angle_param(*elems)  # atom_name_1, atom_name_2, atom_name_3, param_name)
        elif param_type == "dihe":
            if len(elems) != 5:
                quit_wrong_params(line_number, line)
            param = dihe_param(*elems)  # atom_name_1, atom_name_2, atom_name_3, atom_name_4, param_name)
        elif param_type == "branch_dihe":
            if len(elems) != 5:
                quit_wrong_params(line_number, line)
            param = dihe_param(*elems, at_branchpoint=True)
        elif param_type == "branch_imp":
            if len(elems) != 5:
                quit_wrong_params(line_number, line)
            param = dihe_param(*elems, at_branchpoint=True)
        else:  # This can only happen, by programmers error, not by users error
            print("Unknown parameter type: {}.".format(param_type))
        if param.get_param_has_comma():
            self.set_double_params_given(True)
        return param

    def check_entries_lines(self, lines):
        reading_bonds = False
        reading_angles = False
        reading_imp = False
        reading_dihe = False
        reading_branch_diheds = False
        reading_branch_imps = False
        if True in ['ANGLES' in line for line in lines]:
            reading_angles = True
        if True in ['BONDS' in line for line in lines]:
            reading_bonds = True
        if True in ['IMPROPERS' in line for line in lines]:
            reading_imp = True
        if True in ['DIHEDRALS' in line for line in lines]:
            reading_dihe = True
        if True in ['BRANCH DIHEDS' in line for line in lines]:
            reading_branch_diheds = True
            self.set_branch_diheds_given(True)
        if True in ['BRANCH IMPROPERS' in line for line in lines]:
            self.set_branch_imps_given(True)
            reading_branch_imps = True
        return reading_bonds, reading_angles, reading_imp, reading_dihe, reading_branch_diheds, reading_branch_imps

    def read_relevant_lines(self, lines, keyword, param_keyword, list_to_be_appended):
        """This is a template function to read arbitrary sections of the
        listparms file. Crawls through all lines, looks for a specific keyword
        and appends it to specified list."""
        keyword_found = False
        for line_number, line in enumerate(lines, 1):
            if keyword in line:
                keyword_found = True
                continue
            if keyword_found:
                line = strip_comments(line)
                if len(line) == 0:
                    continue
                if line[0].isspace():
                    break  # stop at first blank line
                param = self.create_param(line_number, line, param_keyword)
                list_to_be_appended.append(param)

    def read_bonds_params(self, lines):
        self.read_relevant_lines(lines, keyword="BONDS", param_keyword="bond",
                                 list_to_be_appended=self.__bondparams)

    def read_angle_params(self, lines):
        self.read_relevant_lines(lines, keyword="ANGLES",  param_keyword="angle",
                                 list_to_be_appended=self.__angleparams)

    def read_imp_params(self, lines):
        self.read_relevant_lines(lines, keyword="IMPROPERS",  param_keyword="dihe",
                                 list_to_be_appended=self.__impparams)

    def read_dihe_params(self, lines):
        self.read_relevant_lines(lines, keyword="DIHEDRALS",  param_keyword="dihe",
                                 list_to_be_appended=self.__diheparams)

    def read_branch_dihe_params(self, lines):
        self.read_relevant_lines(lines, keyword="BRANCH DIHEDS",  param_keyword="branch_dihe",
                                 list_to_be_appended=self.__branch_diheparams)

    def read_branch_imp_params(self, lines):
        self.read_relevant_lines(lines, keyword="BRANCH IMPROPERS",  param_keyword="branch_imp",
                                 list_to_be_appended=self.__improper_branch_params)

    def read_bond_params(self, filename):  # This is a very misleading name for this function!
        # read in lines from file
        with open(filename, 'r') as f:
            lines = f.readlines()
        # check for existence of keywords in ffparams file
        switches = self.check_entries_lines(lines)
        reading_bonds, reading_angles, reading_imp, reading_dihe, reading_branch_diheds, reading_branch_imps = switches
        if reading_bonds:
            self.read_bonds_params(lines)
        else:
            print("No bond information in ffparams-file")

        if reading_angles:
            self.read_angle_params(lines)

        if reading_imp:
            self.read_imp_params(lines)

        if reading_dihe:
            self.read_dihe_params(lines)

        if reading_branch_diheds:
            self.read_branch_dihe_params(lines)
        if reading_branch_imps:
            self.read_branch_imp_params(lines)
        if self.get_double_params_given():
            print("!!! WARNING !!!")
            print("You have given multiple parameters for a specific dihedrals.")
            print("Be sure, that you prepared your input properly.")
            print("Defining dihedrals in multiple different ways may lead to more dihedral definitions in your "
                  "topology than you might want to have!\n")

    def print_bond_params(self):
        for param in self.__bondparams:
            print(param.get_atom_1_name() + "   " + param.get_atom_2_name() + "   " + param.get_param())
