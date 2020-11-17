'''
Created on Dec 22, 2013

@author: root
'''
from __future__ import division
from Itp import Itp
from Bond import Bond
import copy
from Utils import ngen
from Utils import core_filename
from Utils import ter_filename
from Utils import inter_filename
from Utils import verbose as Ut_verbose
from BuildingBlockFactory import BuildingBlockFactory
# import time


class Dendrimer(Itp):
    coreFactory = None
    terminalFactory = None
    interFactory = None

    def is_bond(self, atom_1, atom_2):
        for bond in self.get_bond_list():
            if (bond.contains_atom(atom_1) and bond.contains_atom(atom_2) and atom_1.get_nr() != atom_2.get_nr()):
                return True
        return False

    def clear_branches(self):
        self.branch_list = []

    def connectToCore(self, core):
        atoms = copy.copy(core.get_atom_list())
        bonds = copy.copy(core.get_bond_list())
        angles = copy.copy(core.get_angle_list())
        dihedrals = copy.copy(core.get_dihedral_list())
        branches = copy.copy(core.get_branch_list())

        self.set_atom_list(atoms)
        # self.create_coords()
        self.set_bond_list(bonds)
        self.set_angle_list(angles)
        self.set_dihedral_list(dihedrals)
        self.set_branch_list(branches)

    def connectToBb(self, bbFactory):
        branches = copy.copy(self.get_branch_list())
        for atom in self.get_atom_list():
            atom.set_bond_list([])
            atom.set_atom_list([])

        for bond in self.get_bond_list():
            bond.get_a_1().get_atom_list().append(bond.get_other_atom(bond.get_a_1()))
            bond.get_a_2().get_atom_list().append(bond.get_other_atom(bond.get_a_2()))

        for branch_count, branch in enumerate(branches,1):
            if Ut_verbose:
                print("Processing branchpoint nr: {}".format(branch_count))
            bb = bbFactory.makeBBForItp(self)
            # bb is a BuildingBlock NOT a backbone
            donor_branch = bb.get_donor_branch()
            previousList = copy.copy(self.get_atom_list())
            # print(len(previousList))
            self.get_atom_list().extend(bb.get_atom_list())

            acceptor = self.find_atom(branch.get_acceptor())
            # print(type(acceptor), acceptor.get_is_acceptor())  # debug print
            acceptor.set_is_acceptor(True)
            # print(type(acceptor), acceptor.get_is_acceptor())  # debug print
            # print(acceptor.get_nr())
            donor = bb.find_atom(donor_branch.get_donor())
            donor.set_is_donor(True)
            # print(type(donor), donor.get_is_donor(), donor.get_is_branchpoint())  # debug print
            if len(donor.get_atom_list()) < 1 or len(donor.get_atom_list()) > 3:
                print("Your donor atom has a valency of zero or above 2 and we are not prepared to deal with that yet")
                print("Please, contact us, we will be glad to help: pypolybuilder@gmail.com")
                exit()

            # TEST: trying to define bonded terms before:
            bond = Bond(acceptor, donor, 2, None)

            self.get_bond_list().append(bond)
            self.get_branch_list().remove(branch)
            bb.get_branch_list().remove(donor_branch)

            self.get_bond_list().extend(bb.get_bond_list())

            # Useful for v2:
            # The v2 version will have an option to take the angle and dihedral types
            # from the BB file instead of generating them
            self.get_angle_list().extend(bb.get_angle_list())
            self.get_dihedral_list().extend(bb.get_dihedral_list())
            self.get_branch_list().extend(bb.get_branch_list())
            self.get_exclusion_list().extend(bb.get_exclusion_list())

    def read_bbs(self):
        # Inicia topout com os atomos do cor
        # print("Building CORE block")
        self.coreFactory = BuildingBlockFactory(core_filename)
        # print("Building TER block")
        self.terminalFactory = BuildingBlockFactory(ter_filename)
        # print("Building INTER block")
        self.interFactory = BuildingBlockFactory(inter_filename)

        print("Module DENDRIMER was activated")
        print("Reading files:")
        print("Dendrimer core-block: "+ core_filename)
        if not (inter_filename is None):
            print("Dendrimer intr-block: " + inter_filename)
        print("Dendrimer term-block: "+ ter_filename)

    def read_connections(self):
        pass

    def buildAngDihe(self):
        for atom in self.get_atom_list():
            atom.set_bond_list([])
            atom.set_atom_list([])

        for bond in self.get_bond_list():
            bond.get_a_1().get_atom_list().append(bond.get_other_atom(bond.get_a_1()))
            bond.get_a_2().get_atom_list().append(bond.get_other_atom(bond.get_a_2()))

        print("Building Angles...\n")
        self.set_angle_list(self.create_all_angles(self.get_bond_list()))
        print("Building Dihedrals...\n")
        self.set_dihedral_list(self.create_all_dihedrals())

    def connect_bbs(self):
        core = self.coreFactory.makeBBForItp(None)
        self.connectToCore(core)
        #self.buildAngDihe()
        i = 0
        print("Connecting building blocks of generation: " + str(i))

        while (i < ngen):
            i = i + 1
            self.connectToBb(self.interFactory)
            # self.generate_exclusions()
            print("Connecting building blocks of generation: " + str(i))
            #self.buildAngDihe()

        self.connectToBb(self.terminalFactory)
        #self.buildAngDihe()

    def __init__(self):
        self.atom_list = []
        self.bond_list = []
        self.angle_list = []
        self.dihedral_list = []
        self.exclusion_list = []
        self.pair_list = []
        self.set_branch_list([])

        for atom in self.get_atom_list():
            atom.set_bond_list([])
            atom.set_atom_list([])

        for bond in self.get_bond_list():
            bond.get_a_1().get_atom_list().append(bond.get_other_atom(bond.get_a_1()))
            bond.get_a_2().get_atom_list().append(bond.get_other_atom(bond.get_a_2()))

        self.generate_exclusions()
