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
import time
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
    startConnectBB = time.time()
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
    endConnectBB = time.time()
    
    print("Building Angles...")
    startBuildAng = time.time()
    top.set_angle_list(top.create_all_angles(top.get_bond_list()))
    endBuildAng = time.time()
    
    print("Building Dihedrals...")
    startBuildDihed = time.time()
    if Utils.verbose:
        print("Specfific information about specified dihedrals at branchpoints following, eventually: \n")
    # -----------------------  dihedrals  ### ###
    top.set_dihedral_list(top.create_all_dihedrals(forzmatrix=False))
    endBuildDihed = time.time()

    startClear = time.time()
    top.clear_repeated_dihedrals()
    # -----------------------  angles  ### ###
    top.clear_repeated_angles()
    endClear = time.time()
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
    startZMatrix = time.time()
    # angles
    top.set_angle_list(top.create_all_angles(top.get_bond_list()))
    # dihedrals
    top.set_dihedral_list(top.create_all_dihedrals(forzmatrix=True))
    top.clear_repeated_dihedrals()
    top.clear_repeated_angles()
    endZMatrix = time.time()
    
    print("Optimizing Torsions...")
    # quit()
    startOptDihed = time.time()
    optGeo = TorsionOpt(top)
    optGeo.run()
    endOptDihed = time.time()
    
    print("Successfully performed optimization of torsional degrees of freedom!")
    top.calculateAtomsNeighbors()
    top.createPairList()
    top.createNonBondedPairs()
    top.standardPairList()
    
    print("Trying to optimize geometry ...")
    startOptGeom = time.time()
    minGeo = OptimizeGeometry(top)
    minGeo.run()
    endOptGeom = time.time()
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

    if False: # Only used to evaluate performance when needed
        print("***************************")
        print("#------ Time report: ------#")
        print(f'Building angle list: {endBuildAng-startBuildAng:.4f} seconds')
        print(f'Building dihedral list: {endBuildDihed-startBuildDihed:.4f} seconds')
        print(f'Clearing repeated angles and dihedrals: {endClear-startClear:.4f} seconds')
        print(f'Building the Z-Matrix: {endZMatrix-startZMatrix:.4f} seconds')
        print(f'Dihedrals optimization: {endOptDihed-startOptDihed:.4f} seconds')
        print(f'Geometry local optimization: {endOptGeom-startOptGeom:.4f} seconds')
        print("***************************\n")

    print("If your system has any special feature that is not covered by pyPolyBuilder or")
    print("If you encounter any bug, please, help us improve the code at www.github.xxx/pypolybuilder\n")
