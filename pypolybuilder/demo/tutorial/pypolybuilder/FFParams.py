'''
Created on Jan 22, 2015

@author: Vitor
'''
from BondParams  import bond_param
from AngleParams import angle_param
from DihedralParams  import dihe_param

class ff_params(object):
    __bondparams  = None
    __angleparams = None
    __impparams   = None
    __diheparams  = None
    __transdihe   = None
    __cisdihe     = None

    def __init__(self,params_filename):
        self.__bondparams  = []
        self.__angleparams = []
        self.__impparams = []
        self.__diheparams = []
        self.__trans = []
        self.__cis = []
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

    def get_is_trans(self,dihedral):
        a1 = str(dihedral.get_a_1())
        a2 = str(dihedral.get_a_2())
        a3 = str(dihedral.get_a_3())
        a4 = str(dihedral.get_a_4())

        for type in self.__trans:
            if a1 == type[0] and a2 == type[1] and a3 == type[2] and a4 == type[3]:
                return True
        return False

    def get_is_cis(self,dihedral):
        a1 = str(dihedral.get_a_1())
        a2 = str(dihedral.get_a_2())
        a3 = str(dihedral.get_a_3())
        a4 = str(dihedral.get_a_4())

        for type in self.__cis:
            if a1 == type[0] and a2 == type[1] and a3 == type[2] and a4 == type[3]:
                return True
        return False

    def createConformation(self,line):
        line = line.split()
        a1 = line[0]
        a2 = line[1]
        a3 = line[2]
        a4 = line[3]
        conform = line[4]

        if conform == 'trans':
            self.__trans.append([a1,a2,a3,a4])
            return
        self.__cis.append([a1, a2, a3, a4])

    def read_conformations(self,filename):
        f = open(filename, 'r')
        reading_conformation = False
        for line in f:
            if reading_conformation:
                if line.isspace():
                    break
                self.createConformation(line)
            if 'CONFORMATION' in line:
                reading_conformation = True



    def create_param(self,line,param_type):
        atom_name_1 = ''
        atom_name_2 = ''
        atom_name_3 = ''
        atom_name_4 = ''
        param_name       = ''
        i = 0
        
        #reading first atom name
        while(line[i].isspace() == False):
            atom_name_1 = atom_name_1 + line[i]
            i = i + 1
        
        #skipping blank
        while(line[i].isspace()):
            i = i + 1
            
        while(line[i].isspace() == False):
            atom_name_2 = atom_name_2 + line[i]
            i = i + 1
        
        #print(line[i])    
        while(line[i].isspace()):
            i = i + 1
        
        if(param_type == "bond"):
            while(line[i].isspace() == False):
                param_name = param_name + line[i]
                i = i + 1
        
            param = bond_param(atom_name_1,atom_name_2,param_name)
            return param
        
        while(line[i].isspace() == False):
            atom_name_3 = atom_name_3 + line[i]
            i = i + 1
            
        while(line[i].isspace()):
            i = i + 1
        
        if(param_type == "angle"):
            while(line[i].isspace() == False):
                param_name = param_name + line[i]
                i = i + 1
            param = angle_param(atom_name_1,atom_name_2,atom_name_3,param_name)
            return param
        
        while(line[i].isspace() == False):
            atom_name_4 = atom_name_4 + line[i]
            i = i + 1
            
        while(line[i].isspace()):
            i = i + 1
        
        while(line[i].isspace() == False):
            param_name = param_name + line[i]
            i = i + 1
            param = dihe_param(atom_name_1,atom_name_2,atom_name_3,atom_name_4,param_name)
        
        return param
            
    def read_bond_params(self,filename):  
        f = open(filename, 'r')
        param          = None
        reading_bonds  = False
        reading_angles = False
        reading_imp    = False
        reading_dihe   = False
        i = 0
        for line in f:                                              #percorrendo arquivo
            if(reading_bonds and line[0].isspace() == False):       #verificando se esta nos bonds e se nao eh uma linha em branco
                if('ANGLES' in line):                               #verificando se atingimos linha dos angles
                    reading_bonds  = False
                    reading_angles = True
                else:
                    param = self.create_param(line,"bond")
                    self.__bondparams.append(param)
            elif(reading_angles and line[0].isspace() == False):
               if('IMPROPERS' in line):                              
                    reading_angles = False
                    reading_imp    = True
               else:
                    param = self.create_param(line,"angle")
                    self.__angleparams.append(param)
            elif(reading_imp and line[0].isspace() == False):
               if('DIHEDRALS' in line):                              
                    reading_imp  = False
                    reading_dihe = True
               else:
                    param = self.create_param(line,"dihe")
                   # print(param)
                    self.__impparams.append(param)
            elif(reading_dihe and line[0].isspace() == False):
               if('CONFORMATION' in line):
                    reading_dihe = False
               else:
                    param = self.create_param(line,"dihe")
                    self.__diheparams.append(param)
            elif('BONDS' in line):
                reading_bonds = True
                
        #self.print_bond_params()
  
  
    def print_bond_params(self):
        for param in self.__bondparams:
            print(param.get_atom_1_name() + "   " + param.get_atom_2_name() + "   " + param.get_param())
