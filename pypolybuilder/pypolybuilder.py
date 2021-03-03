#!/usr/bin/env python
# !-*- coding: utf8 -*-

from Itp import Itp
import getopt
import sys
import copy
from BuildingBlock import BuildingBlock
from BuildingBlockFactory import BuildingBlockFactory
from Printer import Printer
from Polymer import Polymer
from Dendrimer import Dendrimer
from OptimizeGeometry import OptimizeGeometry
from TorsionOpt import TorsionOpt
import Utils
from constants import DENDRIMER
from constants import POLYMER
from constants import GROMOS_FORMAT
from constants import GROMACS_FORMAT
# import shutil


def main():
    Output = {}
    Output[DENDRIMER] = Dendrimer
    Output[POLYMER] = Polymer
    mode = Utils.buildingMode

    while mode != DENDRIMER and mode != POLYMER:
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
    top.generate_pairs14()
    # catch recursion error of deepcopy
    if len(top.get_atom_list()) > 1000:
        print("You are building a very large structure!\nYour topology will include {} atoms.\n"
              "This will take a while.\n".format(len(top.get_atom_list())))
        sys.setrecursionlimit(100000)
        print("Increasing recursionlimit...\n\n")

    print("Building Angles...")
    top.set_angle_list(top.create_all_angles(top.get_bond_list()))
    print("Building Dihedrals...")
    if Utils.verbose:
        print("Specfific information about specified dihedrals at branchpoints following, eventually: \n")
    # -----------------------  dihedrals  ### ###
    top.set_dihedral_list(top.create_all_dihedrals(forzmatrix=False))
    top.clear_repeated_dihedrals()
    # -----------------------  angles  ### ###
    top.clear_repeated_angles()
    print("Topology is done!\n")
    # -----------------------  write topology ------------------------------ #
    if Utils.outputFormat == GROMOS_FORMAT:
        Printer.printGromos(top)
    elif Utils.outputFormat == GROMACS_FORMAT:
        Printer.printGromacs(top)
    else:
        print("Unkown Output format! Printing Gromacs format instead:")
        Printer.printGromacs(top)

    if Utils.nogeom:
        quit()  # test everything before torsion opt
    # -----------------------  geom opt ------------------------------ #
    print("Building initial structure...")
    # angles
    top.set_angle_list(top.create_all_angles(top.get_bond_list()))
    # dihedrals
    top.set_dihedral_list(top.create_all_dihedrals(forzmatrix=True))
    top.clear_repeated_dihedrals()
    top.clear_repeated_angles()
    print("Optimizing Torsions...")
    # quit()
    optGeo = TorsionOpt(top)
    optGeo.run()
    print("Successfully performed optimization of torsional degrees of freedom!")
    top.calculateAtomsNeighbors()
    top.createPairList()
    top.createNonBondedPairs()
    top.standardPairList()
    print("Trying to optimize geometry ...")
    minGeo = OptimizeGeometry(top)
    minGeo.run()
    # -----------------------  write coordinates ------------------------------ #
    if Utils.outputFormat == GROMOS_FORMAT:
        Printer.printGROMOSCoordFile(top)
    elif Utils.outputFormat == GROMACS_FORMAT:
        Printer.printGROMACSCoordFile(top)
    else:
        print("Unknown Output format! Printing Gromacs format instead:")
        Printer.printGROMACSCoordFile(top)

    print("All Done! Please don't forget to\n")
    print("***************************")
    print("*!!!CHECK YOUR TOPOLOGY!!!*")
    print("***************************\n")

    print("If your system has any special feature that is not covered by pyPolyBuilder or")
    print("If you encounter any bug, please, help us improve the code at www.github.xxx/pypolybuilder\n")
