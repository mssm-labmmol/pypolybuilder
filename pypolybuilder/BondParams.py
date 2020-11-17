'''
Created on Jan 22, 2015

@author: Vitor
'''

class bond_param(object):
    __atom_1_name = None
    __atom_2_name = None
    __param       = None

    def __init__(self, name_1,name_2,param):
        self.__atom_1_name = name_1
        self.__atom_2_name = name_2
        self.__param       = param

    
    
    def get_atom_1_name(self):
        return self.__atom_1_name


    def get_atom_2_name(self):
        return self.__atom_2_name


    def get_param(self):
        return self.__param


    def set_atom_1_name(self, value):
        self.__atom_1_name = value


    def set_atom_2_name(self, value):
        self.__atom_2_name = value


    def set_param(self, value):
        self.__param = value


    def del_atom_1_name(self):
        del self.__atom_1_name


    def del_atom_2_name(self):
        del self.__atom_2_name


    def del_param(self):
        del self.__param

    atom_1_name = property(get_atom_1_name, set_atom_1_name, del_atom_1_name, "atom_1_name's docstring")
    atom_2_name = property(get_atom_2_name, set_atom_2_name, del_atom_2_name, "atom_2_name's docstring")
    param = property(get_param, set_param, del_param, "param's docstring")
    
    
    