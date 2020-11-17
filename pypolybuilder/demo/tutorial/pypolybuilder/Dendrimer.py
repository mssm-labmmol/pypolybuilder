'''
Created on Dec 22, 2013

@author: root
'''
from __future__ import division
from Itp import Itp
from Angle import Angle
from Bond import Bond
from Dihedral import Dihedral
import random
import copy
from Utils import ngen
from Utils import core_filename
from Utils import ter_filename
from Utils import inter_filename
from BuildingBlockFactory import BuildingBlockFactory


class Dendrimer(Itp):
    coreFactory = None
    terminalFacory = None
    interFactory = None

    def is_bond(self, atom_1, atom_2):
        for bond in self.get_bond_list():
            if (bond.contains_atom(atom_1) and bond.contains_atom(atom_2) and atom_1.get_nr() != atom_2.get_nr()):
                return True
        return False

    def clear_branches(self):
        self.branch_list = []

    def create_coords(self):
        #print("* Creating Molecular Topology *")
        for atom in self.get_atom_list():
            x = random.uniform(0, 1)
            y = random.uniform(0, 1)
            z = random.uniform(0, 1)
            atom.set_x(x)
            atom.set_y(y)
            atom.set_z(z)

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
        for branch in branches:
            bb = bbFactory.makeBBForItp(self)
            donor_branch = bb.get_donor_branch()
            self.get_atom_list().extend(bb.get_atom_list())

            acceptor = self.find_atom(branch.get_acceptor())
            donor = bb.find_atom(donor_branch.get_donor())

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
        self.coreFactory = BuildingBlockFactory(core_filename)
        self.terminalFactory = BuildingBlockFactory(ter_filename)
        self.interFactory = BuildingBlockFactory(inter_filename)

        print("Module DENDRIMER was activated")
        print("Reading files:")
        print("Dendrimer core-block: "+ core_filename)
        if inter_filename != None:
            print("Dendrimer intr-block: " + inter_filename)
        print("Dendrimer term-block: "+ ter_filename)


    def read_connections(self):
        pass

    def connect_bbs(self):
        core = self.coreFactory.makeBBForItp(None)
        self.connectToCore(core)

        i = 0

        #print("Connecting Building Blocks")
        while (i < ngen):
            i = i + 1
            self.connectToBb(self.interFactory)

        self.connectToBb(self.terminalFactory)

    def __init__(self):
        self.atom_list = []
        self.bond_list = []
        self.angle_list = []
        self.dihedral_list = []
        self.exclusion_list = []
        self.pair_list = []
        self.set_branch_list([])

        self.create_coords()

        for atom in self.get_atom_list():
            atom.set_bond_list([])
            atom.set_atom_list([])

        for bond in self.get_bond_list():
            bond.get_a_1().get_atom_list().append(bond.get_other_atom(bond.get_a_1()))
            bond.get_a_2().get_atom_list().append(bond.get_other_atom(bond.get_a_2()))

        self.generate_exclusions()
