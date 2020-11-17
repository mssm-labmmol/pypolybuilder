'''
Created on Jan 22, 2015

@author: Vitor
'''


class dihe_param(object):
    __atom_1_name = None
    __atom_2_name = None
    __atom_3_name = None
    __atom_4_name = None
    __param       = None
    __at_branchpoint = None
    __param_has_comma = False
    __unique_identifier = None

    def __init__(self, name_1, name_2, name_3, name_4, param, at_branchpoint=False):
        self.__atom_1_name = name_1
        self.__atom_2_name = name_2
        self.__atom_3_name = name_3
        self.__atom_4_name = name_4
        self.__param       = param
        self.__at_branchpoint = at_branchpoint
        if "," in param:
            self.__param_has_comma = True

    def contains_atom(self, atom_name):
        if(self.atom_1_name == atom_name or self.atom_2_name == atom_name or
                self.atom_3_name == atom_name or self.atom_4_name == atom_name):
            return True

    def get_atom_1_name(self):
        return self.__atom_1_name

    def get_atom_2_name(self):
        return self.__atom_2_name

    def get_atom_3_name(self):
        return self.__atom_3_name

    def get_atom_4_name(self):
        return self.__atom_4_name

    def get_param(self, with_identifier=False):
        if with_identifier:
            return self.__param, self.__unique_identifier
        else:
            return self.__param

    def get_param_has_comma(self):
        return self.__param_has_comma

    def get_unique_identifier(self):
        return self.__unique_identifier

    def set_unique_identifier(self, value):
        self.__unique_identifier = value

    def set_atom_1_name(self, value):
        self.__atom_1_name = value

    def set_atom_2_name(self, value):
        self.__atom_2_name = value

    def set_atom_3_name(self, value):
        self.__atom_3_name = value

    def set_atom_4_name(self, value):
        self.__atom_4_name = value

    def set_param(self, value):
        self.__param = value

    def del_atom_1_name(self):
        del self.__atom_1_name

    def del_atom_2_name(self):
        del self.__atom_2_name

    def del_atom_3_name(self):
        del self.__atom_3_name

    def del_atom_4_name(self):
        del self.__atom_4_name

    def del_param(self):
        del self.__param

    atom_1_name = property(get_atom_1_name, set_atom_1_name, del_atom_1_name, "atom_1_name's docstring")
    atom_2_name = property(get_atom_2_name, set_atom_2_name, del_atom_2_name, "atom_2_name's docstring")
    atom_3_name = property(get_atom_3_name, set_atom_3_name, del_atom_3_name, "atom_3_name's docstring")
    atom_4_name = property(get_atom_4_name, set_atom_4_name, del_atom_4_name, "atom_4_name's docstring")
    param = property(get_param, set_param, del_param, "param's docstring")
