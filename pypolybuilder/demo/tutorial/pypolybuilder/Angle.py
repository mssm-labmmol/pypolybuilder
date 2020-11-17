'''
Created on Dec 22, 2013

@author: root
'''
import Utils

class Angle(object):
    __a1    = None
    __a2    = None
    __a3    = None
    __func  = 2
    __param = None

    def __init__(self, a1, a2, a3, func, param, isDeterminated = False):
        self.__a1 = a1
        self.__a2 = a2
        self.__a3 = a3
        self.__func = func
        self.__param = param
        self.isDeterminated = isDeterminated
        
        #if(self.param == None):
        self.find_angle_param()
        
    @classmethod
    def makeFromBonds(self,b1,b2):
        a1 = b1.get_a_1()
        a2 = b1.get_a_2()
        a3 = b2.get_a_1()
        a4 = b2.get_a_2()
        
            
        if (a1.get_nr() == a3.get_nr()):
            a = Angle(a2, a1, a4,2,None)
        elif(a1.get_nr() == a4.get_nr()):
            a = Angle(a2, a1, a3,2,None)     
        elif(a2.get_nr() == a3.get_nr()):
            a = Angle(a1, a2, a4,2,None) 
        elif(a2.get_nr() == a4.get_nr()):
            a = Angle(a1, a2, a3,2,None) 
        else:
            return None
        
        if(int(a.get_a_1().get_nr()) > int(a.get_a_3().get_nr())):
            aux = a.get_a_3()
            a.set_a_3(a.get_a_1())
            a.set_a_1(aux)
            
        a.find_angle_param()
            
        return a
    
    def find_angle_param(self):
        if (self.get_param() != None):
            return

        for param in Utils.ff_param.get_angle_params():
            if(param.get_atom_2_name() == self.a2.get_atomtype()):
                if(param.get_atom_1_name() == self.a1.get_atomtype() and param.get_atom_3_name() == self.a3.get_atomtype()):
                    self.__param = param.get_param()
                if(param.get_atom_3_name() == self.a1.get_atomtype() and param.get_atom_1_name() == self.a3.get_atomtype()):
                    self.__param = param.get_param()
        
        if(self.get_param() == None):
            print("Warning: Angle " +  str(self.get_a_1().get_atomtype()) + "-" + str(self.get_a_2().get_atomtype()) + "-" +str(self.get_a_3().get_atomtype()) + " has type none")
        
    def get_a_1(self):
        return self.__a1


    def get_a_2(self):
        return self.__a2


    def get_a_3(self):
        return self.__a3


    def get_func(self):
        return self.__func


    def get_param(self):
        return self.__param


    def set_a_1(self, value):
        self.__a1 = value


    def set_a_2(self, value):
        self.__a2 = value


    def set_a_3(self, value):
        self.__a3 = value


    def set_func(self, value):
        self.__func = value


    def set_param(self, value):
        self.__param = value

    def contains_atom(self, atom):
        if (atom.get_nr() == self.get_a_1().get_nr()):
            return True
        elif (atom.get_nr() == self.get_a_2().get_nr()):
            return True
        elif (atom.get_nr() == self.get_a_3().get_nr()):
            return True
        else:
            return False

    a1 = property(get_a_1, set_a_1, None, None)
    a2 = property(get_a_2, set_a_2, None, None)
    a3 = property(get_a_3, set_a_3, None, None)
    func = property(get_func, set_func, None, None)
    param = property(get_param, set_param, None, None)

