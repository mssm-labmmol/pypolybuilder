'''
Created on Dec 22, 2013

@author: root
'''

class Branch(object):
    __donor    = None
    __acceptor = None

    def __init__(self, donor, acceptor):
        self.__donor = donor
        self.__acceptor = acceptor

 
    def get_donor(self):
        return self.__donor


    def get_acceptor(self):
        return self.__acceptor


    def set_donor(self, value):
        self.__donor = value


    def set_acceptor(self, value):
        self.__acceptor = value

    
    donor = property(get_donor, set_donor, None, None)
    acceptor = property(get_acceptor, set_acceptor, None, None)
