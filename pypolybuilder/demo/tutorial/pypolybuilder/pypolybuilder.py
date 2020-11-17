#!/usr/bin/env python
#!-*- coding: utf8 -*-

from Itp       import Itp
import getopt
import sys
import copy
from BuildingBlock import BuildingBlock
from BuildingBlockFactory import BuildingBlockFactory
from Printer import Printer
from Polymer import Polymer
from Dendrimer import Dendrimer
from OptimizeGeometry import OptimizeGeometry
import Utils
from constants import DENDRIMER
from constants import POLYMER
from constants import GROMOS_FORMAT
from constants import GROMACS_FORMAT
import shutil

def main():
    Output = {}
    Output[DENDRIMER] = Dendrimer
    Output[POLYMER] = Polymer
    mode = Utils.buildingMode

    while (mode != DENDRIMER and mode != POLYMER):
        mode = raw_input('Choose 1 for Dendrimer or 2 for Polymer construction: ')

    mode = int(mode)
    top = Output[mode]()

    # Read and create bb_factories dict. Ex: bb_factories[CCO] = factory do CCO
    top.read_bbs()
    # Read and create connections array
    top.read_connections()
    # connect building blocks
    print("Connecting Building Blocks...")
    top.connect_bbs()


    for atom in top.get_atom_list():
        atom.set_bond_list([])
        atom.set_atom_list([])

    for bond in top.get_bond_list():
        bond.get_a_1().get_atom_list().append(bond.get_other_atom(bond.get_a_1()))
        bond.get_a_2().get_atom_list().append(bond.get_other_atom(bond.get_a_2()))


    top.generate_exclusions()

    print("Building Angles...")
    top.set_angle_list(top.create_all_angles(top.get_bond_list()))
    print("Building Dihedrals...\n")
    top.set_dihedral_list(top.create_all_dihedrals(top.get_bond_list()))


    print("All Done! Please don't forget to\n")
    print("***************************")
    print("*!!!CHECK YOUR TOPOLOGY!!!*")
    print("***************************\n")

    print("If your system has any special feature that is not covered by pyPolyBuilder or")
    print("If you encounter any bug, please, help us improve the code at www.github.xxx/pypolybuilder\n")


    # top.calculateAtomsNeighbors()
    #top.create_initial_coords()


    if(Utils.optimizeGeometry):

        #print("Calculating Atom Neighbors")
        top.calculateAtomsNeighbors()
        #print("Creating Pair List")
        top.createPairList()
        #print("Creating Non-bonded pairs")
        top.createNonBondedPairs()
        #print("Standard pair list")
        top.standardPairList()
        print("Trying to optimize geometry ...")
        optGeo = OptimizeGeometry(top)
        optGeo.run()
        print("Done!")
        print("\n***WARNING***: This geometry serves only as an initial guess for further geometry optimization using the MD code of your choice (e.g. GROMOS or GROMACS)")

    top.clear_repeated_dihedrals()
    top.clear_repeated_angles()

    if(Utils.outputFormat == GROMOS_FORMAT):
        Printer.printGromos(top)
        Printer.printGROMOSCoordFile(top)
    else:
        Printer.printGromacs(top)
        Printer.printGROMACSCoordFile(top)
