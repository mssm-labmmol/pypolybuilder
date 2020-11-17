'''
Created on Dec 22, 2013

@author: root
'''

import Utils
import Itp
from Atom import Atom

class Dihedral(object):
    __a1    = None
    __a2    = None
    __a3    = None
    __a4    = None
    __func  = None
    __param = None
    __conformation = None

    def __init__(self, a1, a2, a3,a4, func, param, isDeterminated = False):
        self.__a1 = a1
        self.__a2 = a2
        self.__a3 = a3
        self.__a4 = a4
        self.__func = func
        self.__param = param
        self.isDeterminated = isDeterminated

        if(self.func == "imp"):
            self.find_improper_param()
            self.func = 2
        elif(self.func == "prop"):
            self.find_dihe_param()
            self.func = 1
            if(Utils.ff_param.get_is_trans(self)): self.__conformation = 'trans'
            elif(Utils.ff_param.get_is_cis(self)): self.__conformation = 'cis'

    @classmethod
    def newMake(cls,atom):
        dihe_list = []
        if(len(atom.get_atom_list()) == 3):
            d = Dihedral(atom,atom.get_atom_list()[0],atom.get_atom_list()[1],atom.get_atom_list()[2],"imp",None)
            dihe_list.append(d)
            
        if(len(atom.get_atom_list()) == 2 or len(atom.get_atom_list()) == 3):
            vizinho1 = atom.get_atom_list()[0]
            vizinho2 = atom.get_atom_list()[1]
            
            for vizinhoTemp in vizinho1.get_atom_list():
                if(vizinhoTemp != None and vizinhoTemp != atom):
                    d = Dihedral(vizinhoTemp,vizinho1,atom,vizinho2,"prop",None)
                    dihe_list.append(d)
                    break
                    
            for vizinhoTemp in vizinho2.get_atom_list():
                if(vizinhoTemp != None and vizinhoTemp != atom):
                    d = Dihedral(vizinhoTemp,vizinho2,atom,vizinho1,"prop",None)
                    dihe_list.append(d)
                    break
                
        return dihe_list
        
    def containsAtom(self,atom):
        nr = atom.get_nr()
        if(self.get_a_1().get_nr() == nr or self.get_a_2().get_nr() == nr or self.get_a_3().get_nr() == nr or self.get_a_4().get_nr() == nr):
            return True
        return False
    
    def find_improper_param(self):
        for param in Utils.ff_param.get_improper_params():
            if(param.contains_atom(self.a1.get_atomtype()) and param.contains_atom(self.a2.get_atomtype())
                and param.contains_atom(self.a3.get_atomtype()) and param.contains_atom(self.a4.get_atomtype())):
                     self.__param = param.get_param()

        if(self.get_param() == None):
            print("Warning: Improper Dihedral " +  str(self.get_a_1().get_atomtype()) + "-" + str(self.get_a_2().get_atomtype()) + "-" +str(self.get_a_3().get_atomtype()) +"-"+ str(self.get_a_4().get_atomtype()) + " has type none")

    def find_dihe_param(self):
        for param in Utils.ff_param.get_dihe_params():
            if(param.get_atom_2_name() == self.a2.get_atomtype() and param.get_atom_3_name() == self.a3.get_atomtype()):
                if(param.get_atom_1_name() == self.a1.get_atomtype() and param.get_atom_4_name() == self.a4.get_atomtype()):
                    self.__param = param.get_param()
                elif(param.get_atom_1_name() == self.a4.get_atomtype() and param.get_atom_4_name() == self.a1.get_atomtype()):
                    self.__param = param.get_param()
            if(param.get_atom_2_name() == self.a3.get_atomtype() and param.get_atom_3_name() == self.a2.get_atomtype()):
                if(param.get_atom_1_name() == self.a1.get_atomtype() and param.get_atom_4_name() == self.a4.get_atomtype()):
                    self.__param = param.get_param()
                elif(param.get_atom_1_name() == self.a4.get_atomtype() and param.get_atom_4_name() == self.a1.get_atomtype()):
                    self.__param = param.get_param()

        if(self.get_param() == None):
            print("Warning: Torsional Dihedral " +  str(self.get_a_1().get_atomtype()) + "-" + str(self.get_a_2().get_atomtype()) + "-" +str(self.get_a_3().get_atomtype()) +"-"+ str(self.get_a_4().get_atomtype()) + " has type none")
    
                    
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


    def get_func(self):
        return self.__func


    def get_param(self):
        return self.__param

    def get_conformation(self):
        return self.__conformation

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

    def __str__(self):
        rep = str(self.get_a_1()) + " " + str(self.get_a_2()) + " " + str(self.get_a_3()) + " " + str(self.get_a_4())
        return rep

    a1 = property(get_a_1, set_a_1, None, None)
    a2 = property(get_a_2, set_a_2, None, None)
    a3 = property(get_a_3, set_a_3, None, None)
    a4 = property(get_a_4, set_a_4, None, None)
    func = property(get_func, set_func, None, None)
    param = property(get_param, set_param, None, None)
    conformation = property(get_conformation, set_conformation, None, None)


