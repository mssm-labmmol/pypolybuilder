'''
Created on Dec 22, 2013

@author: root
Edit...
@author: paquoika
'''

import Utils
import Itp
from Atom import Atom
from numpy import arange as np_arange
from Utils import print_pypoly_warning


def check_prop_dihed_atypes(param, atypes, any_atom_is_terminal=False):
    found = False
    if (param.get_atom_1_name() == atypes[0]) and (param.get_atom_2_name() == atypes[1]) \
            and (param.get_atom_3_name() == atypes[2]) and (param.get_atom_4_name() == atypes[3]):
        found = True  # Found forward growth
    if any_atom_is_terminal:  # grow backwards if terminal atom in list
        if (param.get_atom_1_name() == atypes[3]) and (param.get_atom_2_name() == atypes[2]) \
                and (param.get_atom_3_name() == atypes[1]) and (param.get_atom_4_name() == atypes[0]):
            found = True  # Found backward growth
    return found


def check_imp_dihed_atypes(param, atypes,any_atom_is_terminal=False):
    if any_atom_is_terminal:  # for improper dihedrals, it is okay to have terminal atoms
        pass
    return param.contains_atom(atypes[0]) and param.contains_atom(atypes[1])\
           and param.contains_atom(atypes[2]) and param.contains_atom(atypes[3])


def return_params(dihe_key, alist):
    if dihe_key == "imp":
        ff_params_list = Utils.ff_param.get_improper_params()
        check_func = check_imp_dihed_atypes
    elif dihe_key == "prop":
        ff_params_list = Utils.ff_param.get_dihe_params()
        check_func = check_prop_dihed_atypes
    elif dihe_key == "imp_branch":
        ff_params_list = Utils.ff_param.get_improper_branch_params()
        check_func = check_imp_dihed_atypes
    elif dihe_key == "prop_branch":
        ff_params_list = Utils.ff_param.get_branch_diheparams()
        check_func = check_prop_dihed_atypes
    else:
        quit("Unknown dihekey. This is not a users, but a programmers error.")

    out_param = None
    centeratom_name = None
    dihed_identifier = None
    atypes = [a.get_atomtype() for a in alist]
    any_atom_is_terminal = True in [a.is_terminal() for a in alist]
    for param in ff_params_list:
        if check_func(param, atypes, any_atom_is_terminal=any_atom_is_terminal):
            centeratom_name = param.get_atom_3_name()
            out_param, dihed_identifier = param.get_param(with_identifier=True)
    return out_param, centeratom_name, dihed_identifier


def print_special_dihedral_info(params, alist):
    atom_types = [a.get_atomtype() for a in alist]
    print("{} {} {} {}".format(*atom_types))
    print("Dihedral parameters: {}".format(params), "\n")  # debug


def print_verbose_dihedral_message(func, alist, params):
    if func == "imp":
        print("Following improper dihedral was given multiple parameters:")
    elif func == "imp_branch":
        print("Improper dihedral as branchpoint was identified in input:")
    elif func == " prop_branch":
        print("Branchpoint dihedral was identified in input:")
    else:  # "prop" or anything else...
        print("Following dihedral was given multiple parameters:")
    print_special_dihedral_info(params, alist)


def process_impropers(atom, dihe_list, func="imp", verbose=False, forzmatrix=False, dihed_count=None):

    alist = [atom, atom.get_atom_list()[0], atom.get_atom_list()[1], atom.get_atom_list()[2]]
    looked_up_param, centeratom_name, dihed_identifier = return_params(func, alist)
    if looked_up_param is None:  # dihed was not specified
        pass
    else:  # was specified. Now, check whether it was given double
        params = looked_up_param.split(",")
        if len(params)==1:
            d = Dihedral(*alist, func=func, param=params[0], was_input=False, is_multiparam=False,
                         has_siblings=False, dihed_idx=dihed_count)
            dihed_count += 1
            dihe_list.append(d)
        elif len(params) > 1:
            siblings = np_arange(start=dihed_count, stop=dihed_count + len(params))
            if verbose:
                print_verbose_dihedral_message(func, alist, params)
            for param in params:
                d = Dihedral(*alist, func=func, param=param, was_input=True, is_multiparam=True,
                             around_center_atom=atom.get_atomtype() == centeratom_name, siblings=siblings,
                             dihed_idx=dihed_count)
                dihed_count += 1
                dihe_list.append(d)
        if not forzmatrix:
            pass
    return dihe_list, dihed_count


def process_torsionals(atom, dihe_list, func="prop", verbose=False, forzmatrix=False, dihed_count=None):
    vizinho1 = atom.get_atom_list()[0]
    vizinho2 = atom.get_atom_list()[1]
    for vizinhoTemp in vizinho1.get_atom_list():
        alist = [vizinhoTemp, vizinho1, atom, vizinho2]
        if not (vizinhoTemp is None) and (vizinhoTemp != atom):
            looked_up_param, centeratom_name, dihed_identifier = return_params(func, alist)
            if looked_up_param is None:
                if forzmatrix:  # dihed was not specified
                    d = Dihedral(*alist, func, None, dihed_idx=dihed_count)
                    dihed_count += 1
                    dihe_list.append(d)
            else:  # was specified. Now, check whether it was given double
                params = looked_up_param.split(",")
                if len(params) == 1:
                    d = Dihedral(*alist, func=func, param=params[0], was_input=False, is_multiparam=False,
                                 has_siblings=False, dihed_idx=dihed_count)
                    dihed_count += 1
                    dihe_list.append(d)
                    if not forzmatrix:
                        break
                elif len(params) > 1:
                    # dihedrals who belong to the same parameter set:
                    siblings = np_arange(start=dihed_count, stop=dihed_count+len(params))
                    if verbose:
                        print_verbose_dihedral_message(func, alist, params)
                    for param in params:
                        if atom.get_was_central_atom():
                            if not (dihed_identifier in atom.get_was_central_atom_for()):
                                d = Dihedral(*alist, func=func, param=param, was_input=True,
                                             isDetermined=True, is_multiparam=True,
                                             around_center_atom=atom.get_atomtype() == centeratom_name,
                                             has_siblings=True, siblings=siblings, dihed_idx=dihed_count)
                                dihed_count += 1
                                dihe_list.append(d)
                        else:
                            d = Dihedral(*alist, func=func, param=param, was_input=True,
                                         isDetermined=True, is_multiparam=True,
                                         around_center_atom=atom.get_atomtype() == centeratom_name,
                                         has_siblings=True, siblings=siblings, dihed_idx=dihed_count)
                            dihed_count += 1
                            dihe_list.append(d)
                    if not forzmatrix:
                        atom.set_was_central_atom(True)  # was used as central atom for multiparam dihed
                        atom.append_was_central_atom_for(dihed_identifier)
                        break

    for vizinhoTemp in vizinho2.get_atom_list():
        alist = [vizinhoTemp, vizinho2, atom, vizinho1]
        if not (vizinhoTemp is None) and vizinhoTemp != atom:
            looked_up_param, centeratom_name, dihed_identifier = return_params(func, alist)
            if looked_up_param is None:
                if forzmatrix:  # dihed was not specified
                    d = Dihedral(*alist, func, None, dihed_idx=dihed_count)
                    dihed_count += 1
                    dihe_list.append(d)
            else:  # was specified. Now, check whether it was given double
                params = looked_up_param.split(",")
                if len(params) == 1:
                    d = Dihedral(*alist, func=func, param=params[0], was_input=False, is_multiparam=False,
                                 has_siblings=False, dihed_idx=dihed_count)
                    dihed_count += 1
                    dihe_list.append(d)
                    if not forzmatrix:
                        break
                elif len(params) > 1:
                    # dihedrals who belong to the same parameter set:
                    siblings = np_arange(start=dihed_count, stop=dihed_count+len(params))
                    if verbose:
                        print_verbose_dihedral_message(func, alist, params)
                    for param in params:
                        if atom.get_was_central_atom():
                            if not (dihed_identifier in atom.get_was_central_atom_for()):
                                d = Dihedral(*alist, func=func, param=param, was_input=True,
                                             isDetermined=True, is_multiparam=True,
                                             around_center_atom=atom.get_atomtype() == centeratom_name,
                                             has_siblings=True, siblings=siblings, dihed_idx=dihed_count)
                                dihed_count += 1
                                dihe_list.append(d)
                        else:
                            d = Dihedral(*alist, func=func, param=param, was_input=True,
                                         isDetermined=False, is_multiparam=True,
                                         around_center_atom=atom.get_atomtype() == centeratom_name,
                                         has_siblings=True, siblings=siblings, dihed_idx=dihed_count)
                            dihed_count += 1
                            dihe_list.append(d)
                    if not forzmatrix:
                        atom.set_was_central_atom(True)  # was used as central atom for multiparam dihed
                        atom.append_was_central_atom_for(dihed_identifier)
                        break
    return dihe_list, dihed_count


def process_tetrahedral_center(atom, dihe_list, verbose=False, forzmatrix=False):
    """This function is a legacy, since it is used nowhere!"""
    pass
    return dihe_list


def process_branchpoint_atom(atom, dihe_list, forzmatrix=False, dihed_count=None):
    """This is a very fancy function, which may not be used very often..."""
    alist = atom.get_atom_list()
    if len(alist) == 1:
        # THIS SHOULD ACTUALLY NEVER HAPPEN FOR A BRANCH POINT
        pass

    if len(alist) == 3:
        dihe_list, dihed_count = process_impropers(atom, dihe_list, func="imp_branch",
                                                   verbose=Utils.verbose, forzmatrix=forzmatrix,
                                                   dihed_count=dihed_count)

    if len(atom.get_atom_list()) == 2 or len(atom.get_atom_list()) == 3:
        dihe_list, dihed_count = process_torsionals(atom, dihe_list, func="prop_branch",
                                                    verbose=Utils.verbose, forzmatrix=forzmatrix,
                                                    dihed_count=dihed_count)
    if len(alist) == 4:  # tetrahedral center
        pass
    return dihe_list, dihed_count


def process_non_branchpoint_atom(atom, dihe_list, forzmatrix=False, dihed_count=None):
    """This function is almost the same as process_branchpoint_atom. But only almost..."""
    alist = atom.get_atom_list()
    if len(alist) == 1:  # terminal atom
        pass

    if len(alist) == 3:
        dihe_list, dihed_count = process_impropers(atom, dihe_list, func="imp",
                                                   verbose=Utils.verbose, forzmatrix=forzmatrix,
                                                   dihed_count=dihed_count)

    if len(atom.get_atom_list()) == 2 or len(atom.get_atom_list()) == 3:
        dihe_list, dihed_count = process_torsionals(atom, dihe_list, func="prop",
                                                    verbose=Utils.verbose, forzmatrix=forzmatrix,
                                                    dihed_count=dihed_count)

    if len(alist) == 4:  # tetrahedral center
        pass
        # dihe_list = process_tetrahedral_center(atom, dihe_list, verbose=Utils.verbose, forzmatrix=forzmatrix)

    return dihe_list, dihed_count


# INTRODUCE CLASS #########################################################
class Dihedral(object):
    __a1    = None
    __a2    = None
    __a3    = None
    __a4    = None
    __func  = None
    __param = None
    __conformation = None
    __was_input = None
    __is_multiparam = None
    __around_center_atom = None
    __dihed_idx = None
    __has_siblings = None
    __siblings = []

    def __init__(self, a1, a2, a3, a4, func, param, isDetermined=False, was_input=False,
                 around_center_atom=False, is_multiparam=False, dihed_idx=None, has_siblings=False,
                 siblings=[]):
        self.__a1 = a1
        self.__a2 = a2
        self.__a3 = a3
        self.__a4 = a4
        self.__func = func
        self.__param = param
        self.isDetermined = isDetermined
        self.__was_input = was_input
        self.__a1.connected_list.add(self.__a2)
        self.__a2.connected_list.add(self.__a1)
        self.__a2.connected_list.add(self.__a3)
        self.__a3.connected_list.add(self.__a2)
        self.__a1.connected_list.add(self.__a3)
        self.__a3.connected_list.add(self.__a1)
        self.__a4.connected_list.add(self.__a3)
        self.__a3.connected_list.add(self.__a4)
        self.__a1.connected_list.add(self.__a4)
        self.__a4.connected_list.add(self.__a1)
        self.__a2.connected_list.add(self.__a4)
        self.__a4.connected_list.add(self.__a2)
        self.__is_multiparam = is_multiparam
        self.__around_center_atom = around_center_atom
        self.__dihed_idx = dihed_idx
        self.__has_siblings = has_siblings
        self.__siblings = siblings
        # this asserts are very unusual... I dont think this would ever be found in an itp file...
        if self.func == "imp":
            if self.__param is None:
                self.find_improper_param()
            self.func = 2

        elif self.func == "prop":
            if self.__param is None:
                self.find_dihe_param()
            self.func = 1
            if Utils.ff_param.get_is_trans(self):
                self.__conformation = 'trans'
            elif Utils.ff_param.get_is_cis(self):
                self.__conformation = 'cis'

        elif self.func == "prop_branch":
            if self.__param is None:
                self.find_branch_dihe_param()
            self.func = 1
            if Utils.ff_param.get_is_trans(self):
                self.__conformation = 'trans'
            elif Utils.ff_param.get_is_cis(self):
                self.__conformation = 'cis'

        elif self.func == "imp_branch":
            if self.__param is None:
                self.find_improper_branch_param()
            self.func = 2

    @classmethod
    def newMake(cls, atom, forzmatrix=False, dihed_count=None):
        dihe_list = []
        # Case: dihedral was specified in connect-file. (could be multiple values for same atomset)
        if atom.get_is_in_specified_dihedral():
            if Utils.verbose:
                print("Atom in specified dihedral: {}".format(atom.get_nr()))
            dihed_vals = atom.get_specified_dihedral_vals()
            if dihed_vals["use_for_lookup"]:
                for funct, ph0 in zip(dihed_vals["funct"], dihed_vals["ph0"]):
                    dihe_list.append(Dihedral(*dihed_vals["dihed_atoms"], func=funct, param=ph0, was_input=True,
                                              dihed_idx=dihed_count))
                    dihed_count += 1
        # Case: dihedral was specified in list-parm file. (could be multiple values for same atomset)
        if atom.get_is_branchpoint():
            if Utils.verbose and not forzmatrix:
                print("Atom is at branchpoint: {}".format(atom.get_nr()))
            if Utils.ff_param.get_branch_diheds_given() or Utils.ff_param.get_branch_imps_given():
                # this function breaks the loop eventually:
                dihe_list, dihed_count = process_branchpoint_atom(atom, dihe_list, forzmatrix=forzmatrix,
                                                                  dihed_count=dihed_count)  # static function 
        # either break before, or continue here:
        dihe_list, dihed_count = process_non_branchpoint_atom(atom, dihe_list, forzmatrix=forzmatrix,
                                                              dihed_count=dihed_count)
                
        return dihe_list, dihed_count

    def containsAtom(self, atom):
        nr = atom.get_nr()
        if(self.get_a_1().get_nr() == nr or self.get_a_2().get_nr() == nr or
           self.get_a_3().get_nr() == nr or self.get_a_4().get_nr() == nr):
            return True
        return False

    def isAroundPair(self, atom1, atom2):
        nra1 = atom1.get_nr()
        nra2 = atom2.get_nr()
        if((self.get_a_2.get_nr() == nra1 and self.get_a_3 == nra2) or
           (self.get_a_2.get_nr() == nra2 and self.get_a_3 == nra1)):
            return True
        return False
    
    def find_improper_param(self):
        for param in Utils.ff_param.get_improper_params():
            if(param.contains_atom(self.a1.get_atomtype()) and param.contains_atom(self.a2.get_atomtype())
                    and param.contains_atom(self.a3.get_atomtype()) and param.contains_atom(self.a4.get_atomtype())):
                self.__param = param.get_param()

        if self.get_param() is None:
            atom_types = [a.get_atomtype() for a in [self.get_a_1(), self.get_a_2(), self.get_a_3(), self.get_a_4()]]
            print_pypoly_warning("Warning: Improper Dihedral {}-{}-{}-{} has type none.".format(*atom_types))

    def find_dihe_param(self):
        for param in Utils.ff_param.get_dihe_params():
            if(param.get_atom_2_name() == self.a2.get_atomtype() and
               param.get_atom_3_name() == self.a3.get_atomtype()):
                if(param.get_atom_1_name() == self.a1.get_atomtype() and
                   param.get_atom_4_name() == self.a4.get_atomtype()):
                    self.__param = param.get_param()
                elif(param.get_atom_1_name() == self.a4.get_atomtype() and
                     param.get_atom_4_name() == self.a1.get_atomtype()):
                    self.__param = param.get_param()

            if(param.get_atom_2_name() == self.a3.get_atomtype() and
               param.get_atom_3_name() == self.a2.get_atomtype()):
                if(param.get_atom_1_name() == self.a1.get_atomtype() and
                   param.get_atom_4_name() == self.a4.get_atomtype()):
                    self.__param = param.get_param()
                elif(param.get_atom_1_name() == self.a4.get_atomtype() and
                     param.get_atom_4_name() == self.a1.get_atomtype()):
                    self.__param = param.get_param()

        if self.get_param() is None:
            atom_types = [a.get_atomtype() for a in [self.get_a_1(), self.get_a_2(), self.get_a_3(), self.get_a_4()]]
            print_pypoly_warning("Warning: Torsional Dihedral {}-{}-{}-{} has type none.".format(*atom_types))

    def find_branch_dihe_param(self):
        success = False  # debug
        for param in Utils.ff_param.get_branch_diheparams():
            if(param.get_atom_2_name() == self.a2.get_atomtype() and
                    param.get_atom_3_name() == self.a3.get_atomtype()):
                if(param.get_atom_1_name() == self.a1.get_atomtype() and
                   param.get_atom_4_name() == self.a4.get_atomtype()):
                    self.__param = param.get_param()
                    success = True  # debug
                elif(param.get_atom_1_name() == self.a4.get_atomtype() and
                     param.get_atom_4_name() == self.a1.get_atomtype()):
                    self.__param = param.get_param()
                    success = True  # debug
            if(param.get_atom_2_name() == self.a3.get_atomtype() and
               param.get_atom_3_name() == self.a2.get_atomtype()):
                if(param.get_atom_1_name() == self.a1.get_atomtype() and
                   param.get_atom_4_name() == self.a4.get_atomtype()):
                    self.__param = param.get_param()
                    success = True  # debug
                elif(param.get_atom_1_name() == self.a4.get_atomtype() and
                     param.get_atom_4_name() == self.a1.get_atomtype()):
                    self.__param = param.get_param()
                    success = True  # debug

    def find_improper_branch_param(self):
        for param in Utils.ff_param.get_improper_branch_params():
            if(param.contains_atom(self.a1.get_atomtype()) and param.contains_atom(self.a2.get_atomtype())
                    and param.contains_atom(self.a3.get_atomtype()) and param.contains_atom(self.a4.get_atomtype())):
                self.__param = param.get_param()

        if self.get_param() is None:
            atom_types = [a.get_atomtype() for a in [self.get_a_1(), self.get_a_2(), self.get_a_3(), self.get_a_4()]]
            print_pypoly_warning("Warning: Improper Dihedral {}-{}-{}-{} has type none.".format(*atom_types))

    def print_atoms(self):
        print(str(self.a1.get_nr()) + " " + str(self.a2.get_nr()) + " " + str(self.a3.get_nr()) + " " + str(self.a4.get_nr()) + " ")
        
    def get_a_1(self):
        return self.__a1
    
    def get_a_2(self):
        return self.__a2

    def get_a_3(self):
        return self.__a3
    
    def get_a_4(self):
        return self.__a4

    def get_all_atoms(self):
        return [self.get_a_1(), self.get_a_2(), self.get_a_3(), self.get_a_4()]

    def get_all_atoms_nr(self):
        return [self.get_a_1().get_nr(), self.get_a_2().get_nr(), self.get_a_3().get_nr(), self.get_a_4().get_nr()]

    def is_torsional(self):
        return self.get_func() == 1

    def get_func(self):
        return self.__func

    def get_param(self):
        return self.__param

    def get_conformation(self):
        return self.__conformation

    def get_was_input(self):
        return self.__was_input

    def get_is_multiparam(self):
        return self.__is_multiparam

    def get_around_center_atom(self):
        return self.__around_center_atom

    def get_has_siblings(self):
        return self.__has_siblings

    def get_siblings(self):
        return self.__siblings

    def get_dihed_idx(self):
        return self.__dihed_idx

    def set_is_multiparam(self, value):
        self.__is_multiparam = value

    def set_conformation(self, value):
        self.__conformation = value

    def set_a_1(self, value):
        self.__a1 = value

    def set_a_2(self, value):
        self.__a2 = value

    def set_a_3(self, value):
        self.__a3 = value
    
    def set_a_4(self, value):
        self.__a4 = value

    def set_func(self, value):
        self.__func = value

    def set_param(self, value):
        self.__param = value

    def set_was_input(self, value):
        self.__was_input = value

    def __str__(self):
        alist = [self.get_a_1(), self.get_a_2(), self.get_a_3(), self.get_a_4()]
        rep = alist[0].get_atom() + str(alist[0].get_nr()) + " " \
              + alist[1].get_atom() + str(alist[1].get_nr()) + " " \
              + alist[2].get_atom() + str(alist[2].get_nr()) + " "\
              + alist[3].get_atom() + str(alist[3].get_nr()) + " - " + str(self.get_func())
        return rep

    a1 = property(get_a_1, set_a_1, None, None)
    a2 = property(get_a_2, set_a_2, None, None)
    a3 = property(get_a_3, set_a_3, None, None)
    a4 = property(get_a_4, set_a_4, None, None)
    func = property(get_func, set_func, None, None)
    param = property(get_param, set_param, None, None)
    conformation = property(get_conformation, set_conformation, None, None)
    was_input = property(get_was_input, set_was_input, None, None)
    # is_multiparam = property(get_was_input, set_was_input, None, None)


