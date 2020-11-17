'''
Created on Dec 22, 2013

@author: root
'''
import Utils

class Bond(object):
    __a1    = None
    __a2    = None
    __func  = None
    __param = None
    
    def __init__(self, a1, a2, func, param):
        self.__a1 = a1
        self.__a2 = a2
        self.__func = func
        self.__param = param
        
        #if(self.param == None):
        self.find_bond_param()
    
   
    
    def find_bond_param(self):
        for param in Utils.ff_param.get_bond_params():
            if(param.get_atom_1_name() == self.a1.get_atomtype() and param.get_atom_2_name() == self.a2.get_atomtype()):
                self.__param = param.get_param()
            elif(param.get_atom_2_name() == self.a1.get_atomtype() and param.get_atom_1_name() == self.a2.get_atomtype()):
                self.__param = param.get_param()

        if(self.get_param() == None):  # should this be in the loop?
            print("Warning: Bond " +  str(self.get_a_1().get_atomtype()) + "-" + str(self.get_a_2().get_atomtype()) +  " has type none")

    def get_a_1(self):
        return self.__a1

    def get_a_2(self):
            return self.__a2

    def get_func(self):
            return self.__func

    def get_param(self):
        return self.__param

    def set_a_1(self, value):
        self.__a1 = value

    def set_a_2(self, value):
        self.__a2 = value

    def set_func(self, value):
            self.__func = value

    def set_param(self, value):
            self.__param = value

    def contains_atom(self, atom):
        if(atom.get_nr() == self.get_a_1().get_nr()):
            return True
        elif(atom.get_nr() == self.get_a_2().get_nr()):
            return True
        else:
            return False

    def get_other_atom(self,atom):
        if(atom.get_nr() == self.get_a_1().get_nr()):
            return self.get_a_2()
        else:
            return self.get_a_1()

    def print_bond(self):
        print("a1: " + str(self.a1.get_nr()) + " a2: " + str(self.a2.get_nr()))

        a1 = property(get_a_1, set_a_1, None, None)
        a2 = property(get_a_2, set_a_2, None, None)
        func = property(get_func, set_func, None, None)
        param = property(get_param, set_param, None, None)
