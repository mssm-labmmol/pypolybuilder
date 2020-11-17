'''
Created on Apr 2, 2015

@author: vitor
'''
from Itp import Itp
from Atom     import Atom
from Bond     import Bond
from Angle    import Angle
from Dihedral import Dihedral
from Branch   import Branch
from BuildingBlock import BuildingBlock
import copy

class BuildingBlockFactory(object):

    __template = None


    #Initializes the BuildingBlock setting his template
    def __init__(self, path):
        self.set_template(BuildingBlock(path))
        self.resName = self.template.resName

    def get_template(self):
        return self.__template


    def set_template(self, value):
        self.__template = value


    def del_template(self):
        del self.__template


    #Create a buildingblock for an Itp
    def makeBBForItp(self,itp):
        buildingBlock = copy.deepcopy(self.__template)

        if(itp != None):
            length = len(itp.get_atom_list())
            maxCgnr = 0
            if(length > 0):
                maxCgnr = max(int(node.get_cgnr())for node in itp.get_atom_list())
            for branch in buildingBlock.get_branch_list():
                if(branch.get_donor() != None):
                    branch.set_donor(int(branch.get_donor()) + length) 
                if(branch.get_acceptor() != None):
                    branch.set_acceptor(int(branch.get_acceptor()) + length) 

            for atom in buildingBlock.get_atom_list():
                atom.set_nr(int(atom.get_nr()) + length)
                atom.set_cgnr(int(atom.get_cgnr()) + int(maxCgnr))

            
        return buildingBlock
            
    template = property(get_template, set_template, del_template, "template's docstring")
        
