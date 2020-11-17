"""
Created on Jan 23, 2015

@author: Vitor

Edit...
@author: paquoika
"""

import numpy as np
from pypoly_io import parse_arguments, process_arguments, disentangle_args, debug_print_args
import warnings
import sys

global DEBUG
DEBUG=False


# Throw and format warnings:
def get_filename_for_warnings(stacklevel=1):
    """This is a monkey patch from the warnings module
    https://hg.python.org/cpython/file/2.7/Lib/warnings.py#l37"""
    # Get context information
    try:
        caller = sys._getframe(stacklevel)
    except ValueError:
        globals = sys.__dict__
        lineno = 1
    else:
        globals = caller.f_globals
        lineno = caller.f_lineno
    filename = globals.get('__file__')
    if filename:
        fnl = filename.lower()
        if fnl.endswith((".pyc", ".pyo")):
            filename = filename[:-1]
    return filename, lineno


def print_pypoly_warning(warn_msg, category=RuntimeWarning):
    filename, lineno = get_filename_for_warnings()
    fwarn = warnings.formatwarning(warn_msg, category=category, filename=filename, lineno=lineno, )
    fwarn = "\n".join(fwarn.split("\n")[:-2])  # does not have siblings
    print(fwarn, "\n", file=sys.stderr)


def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def vectorFromAtoms(atom1, atom2):
    vec12 = [atom1.get_x() - atom2.get_x(), atom1.get_y() - atom2.get_y(),
             atom1.get_z() - atom2.get_z()]
    return vec12


def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def deg_to_rad(ang_rad):
    return ang_rad * np.pi / 180.0


def get_ff_bonds(forcefield_path):
    ff_bonds={}
    ffbonded_file = open(forcefield_path+'/ffbonded.itp','r')
    for line in ffbonded_file:
        #Definition of params in gromacs format follows: '#define gb_N B0[N] K[N]'
        if all( [check in line for check in ["#define", "gb_"]] ):
            bond_type = line.split()[1]
            if bond_type not in ff_bonds:
                b0 = float(line.split()[2])
                k  = float(line.split()[3])/1000
                ff_bonds[bond_type]=(b0, k)
            else:
                #This should never happen
                print(f'The bond type {bond_type} is a duplicate.')
    ffbonded_file.close()

    return ff_bonds


def get_ff_angles(forcefield_path):
    ff_angles={}
    ffbonded_file = open(forcefield_path+'/ffbonded.itp','r')
    for line in ffbonded_file:
        #Definition of params in gromacs format follows: '#define ga_N T0[N] K[N]'
        if all( [check in line for check in ["#define", "ga_"]] ):
            angle_type = line.split()[1]
            if angle_type not in ff_angles:
                t0 = float(line.split()[2])
                k  = float(line.split()[3])*10
                ff_angles[angle_type]=(t0, k)
            else:
                #This should never happen
                print(f'The bond type {angle_type} is a duplicate.')
    ffbonded_file.close()

    return ff_angles


def get_ff_dihedrals(forcefield_path):
    ff_dihedrals={}
    ffbonded_file = open(forcefield_path+'/ffbonded.itp','r')
    for line in ffbonded_file:
        #Definition of params in gromacs format follows: '#define gd_N Pd[N] CP[N] NP[N]'
        if all( [check in line for check in ["#define", "gd_"]] ):
            dihed_type = line.split()[1]
            if dihed_type not in ff_dihedrals:
                pd = float(line.split()[2])
                cp = float(line.split()[3])
                np = float(line.split()[4])
                ff_dihedrals[dihed_type]=(pd, cp, np)
            else:
                #This should never happen
                print(f'The bond type {dihed_type} is a duplicate.')
    ffbonded_file.close()

    return ff_dihedrals


def get_ff_LJs(forcefield_path):
    ff_LJs={}
    ffnonbonded_file = open(forcefield_path+'/ffnonbonded.itp','r')
    pairTypes=False
    for line in ffnonbonded_file:
        #Definition of params in gromacs format follows: 'atom1 atom2 func C6 C12'
        if pairTypes and line.split()[0] != ";":
            atom1=line.split()[0]
            atom2=line.split()[1]
            bond_c6=float(line.split()[3])
            bond_c12=float(line.split()[4])*10
            ff_LJs[atom1+" "+atom2]=(bond_c6, bond_c12)
        elif 'pairtypes' in line:
            pairTypes=True
    ffnonbonded_file.close()

    return ff_LJs



def get_ff_params(forcefield_path):
    if forcefield_path == None:
        return None, None, None, None
    else:
        ff_bonds     = get_ff_bonds(forcefield_path)
        ff_angles    = get_ff_angles(forcefield_path)
        ff_dihedrals = get_ff_dihedrals(forcefield_path)
        ff_LJs       = get_ff_LJs(forcefield_path)

    return ff_bonds, ff_angles, ff_dihedrals, ff_LJs


# DIRTY HACK to outsource io, without changing other modules:
# This should better be removed from here, and parse_arguments should be called in __main__
# But this would lead to the necessity to restructure many classes (eg Printer...)
args = parse_arguments()
ff_param, buildingMode, outputFormat = process_arguments(args)

# disentangle:
[output_filename, core_filename, ter_filename, inter_filename,
 gro_filename, ngenga, ngen, nsteps, npop, connections_paths,
 topName, bbs_paths, debug, params_filename, verbose, nogeom,
 forcefield_path, nskipLJ, stepLength] = disentangle_args(args)
if nogeom:
    print("You set the nogeom-flag. Will quit after the building of the topology.\n")
if args.debug:
    debug_print_args(args)

    DEBUG=True

#################################################
