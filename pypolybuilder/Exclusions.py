'''
Created on Apr 20, 2015

@author: vitor
'''

class Exclusions(object):
    '''
    classdocs
    '''
    __exclusion_list = None
    __exclusion_extra = None
    

    def __init__(self, exclusion_list, exclusion_extra):
        exclusion_list.sort()
        self.__exclusion_list = exclusion_list
        self.__exclusion_extra = exclusion_extra

    def print_exclusions(self):
        for exc in self.get_exclusion_list():
            print(exc)
    
    def get_exclusion_list(self):
        return self.__exclusion_list
    
    def get_exclusion_extra(self): #This is used for the Gromacs output only since GROMOS will handle explicitly all the exclusions
        return self.__exclusion_extra

    def set_exclusion_list(self, value):
        self.__exclusion_list = value


    def del_exclusion_list(self):
        del self.__exclusion_list

    exclusion_list = property(get_exclusion_list, set_exclusion_list, del_exclusion_list, "exclusion_list's docstring")

    