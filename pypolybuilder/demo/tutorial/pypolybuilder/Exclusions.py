'''
Created on Apr 20, 2015

@author: vitor
'''

class Exclusions(object):
    '''
    classdocs
    '''
    __exclusion_list = None

    def __init__(self, exclusion_list):
        exclusion_list.sort()
        self.__exclusion_list = exclusion_list

    def print_exclusions(self):
        for exc in self.get_exclusion_list():
            print(exc)
    
    def get_exclusion_list(self):
        return self.__exclusion_list

    def set_exclusion_list(self, value):
        self.__exclusion_list = value


    def del_exclusion_list(self):
        del self.__exclusion_list

    exclusion_list = property(get_exclusion_list, set_exclusion_list, del_exclusion_list, "exclusion_list's docstring")

    