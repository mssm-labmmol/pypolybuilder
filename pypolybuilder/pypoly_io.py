"""
Created on May, 2020

@author: paquoika

"""

import argparse
import os.path
import sys
import shutil
from constants import DENDRIMER
from constants import POLYMER
from constants import GROMOS_FORMAT
from constants import GROMACS_FORMAT
from FFParams import ff_params
from FFParams import give_unique_identifier_to_dihedral_params


def quitError(message):
    print("An error occurred.")
    print(message)
    print("See pypolybuilder -h for more instructions.\n")
    sys.tracebacklimit = 0  # dont show the traceback of the error, when exiting like this
    raise RuntimeError("Initialization failed...")


def check_file_exists(fname):
    if not os.path.exists(fname):
        quitError("{} does not exist.".format(fname))
    else:
        return True


def print_logo():
    script_dir = os.path.dirname(__file__)

    logo_filename = script_dir + '/asciilogo.txt'
    with open(logo_filename, "r") as f:
        shutil.copyfileobj(f, sys.stdout)


def parse_arguments():
    print_logo()

    parser = argparse.ArgumentParser(description='All available options are shown below')
    parser.add_argument('-c', '--core',
                        help='Core filename', dest="core_filename", metavar="file_name.itp")
    parser.add_argument('-i', '--inter',
                        help='Inter filename', default=None, dest="inter_filename", metavar="file_name.itp")
    parser.add_argument('-t', '--ter',
                        help='Ter filename', dest="ter_filename", metavar="file_name.itp",
                        default=None, type=str)
    parser.add_argument('-o', '--output',
                        help='Output filename. Default: default.itp', metavar="file_name.itp",
                        default='default.itp', type=str, dest="output_filename")
    parser.add_argument('-l', '--params',
                        help='Forcefield parameters filename.',
                        dest="params_filename", metavar="file_name.itp", type=str)
    parser.add_argument('-n', '--ngen',
                        help='Number of generations to be built.', default=0, type=int, metavar="ngen")
    parser.add_argument('-nsteps', '--nsteps',
                        help='Number of optimize geometry steps', default=200, type=int, metavar="nsteps")
    parser.add_argument('-ngenga', '--ngenga',
                        help='Number of GA generations', default=20, type=int, metavar="ngenga")
    parser.add_argument('-npop', '--npop',
                        help='Number of GA population.', default=25, type=int, metavar="npop")
    parser.add_argument('-ff', '--forcefield', default=None, type=str, dest='forcefield_path',
                        metavar='path/to/the/force/field', help="The path to get the force field parameter values")
    parser.add_argument('-nskipLJ', '--nskipLJ', default=0, type=int, dest='nskipLJ', metavar="nskipLJ", 
                        help="Number of minimization steps that will be performed without evaluating Lennard-Jones interactions")
    parser.add_argument('-stepLength', '--stepLength', default=0.0001, type=float, dest='stepLength', metavar="Step length", 
                        help="Length of the step used in the geometry optimization")
    # parser.add_argument('-optgeo', '--optgeo',
    #                     help='If enabled pyPolyBuilder will try to find a reasonable geometry. '
    #                          'Switched off by default.',
    #                     dest="optimizeGeometry", action='store_true')

    bmode_args = parser.add_argument_group(title="Building Mode",
                                           description="In which building mode shall pyPolyBuilder operate?")
    bmode_args.add_argument('-dendrimer', '--dendrimer', action='store_true',
                            help="Set building mode to Dendrimer. (This is default.)")
    bmode_args.add_argument('-polymer', '--polymer', action='store_true',
                            help="Set building mode to Polymer.")

    # parser.add_argument('-bbs', '--bbs',
    #                     help='Building blocks filename used to build a Polymer',
    #                     nargs="*", type=list)
    parser.add_argument('-bbs', '--bbs',
                        help='Building blocks filename used to build a Polymer. Provide comma-separated list.',
                        dest="bbs_paths", type=(lambda s: [item for item in s.split(',')]), default=[],
                        metavar="fname1.itp,fname2.itp,...")
    parser.add_argument('-in', '--in',
                        help='Building blocks connections filename', dest="connections_paths",
                        default=None, metavar="file_name.in")
    parser.add_argument('-g', '--gro',
                        help='Gro filename.', default='default.gro',
                        type=str, dest="gro_filename", metavar="file_name.gro")

    outformat_args = parser.add_argument_group(title="Output Format",
                                               description="In which format shall pyPolyBuilder output?")
    outformat_args.add_argument('-gromacs', '--gromacs',
                                help='Set output format to GROMACS. (This is default.)', action='store_true')
    outformat_args.add_argument('-gromos', '--gromos',
                                help='Set output format to GROMOS.', action='store_true')

    parser.add_argument('-name', '--name',
                        help='Topology name to be printed on output file.', dest="topName", default=None,
                        metavar="file_name")
    parser.add_argument('-v', '--verbose',
                        help='Be verbose.', action="store_true")
    parser.add_argument('--debug',
                        help='Print debug messages.', action="store_true")
    parser.add_argument("--nogeom",
                        help="do not perform geometry optimization",
                        action="store_true")

    args = parser.parse_args()
    sanity_check_args(args)

    return args


def sanity_check_buildingode(building_mode, args):
    if building_mode == POLYMER:
        if len(args.bbs_paths) < 1:
            quitError("Please specify an existing polymer building blocks by using --bbs option.")
        if args.connections_paths is None:
            quitError("Please specify a polymer connections filename by using --in option.")
        check_file_exists(args.connections_paths)

    if building_mode == DENDRIMER:
        if args.core_filename is None:
            quitError("Please specify an existing dendrimer core filename by using --core option.")
        check_file_exists(args.core_filename)
        if args.ter_filename is None:
            quitError("Please specify an existing dendrimer terminal filename by using --ter option.")
        check_file_exists(args.ter_filename)


def sanity_check_args(args):

    if args.params_filename is None:
        quitError("Please specify an existing Force Field parameters filename by using --params option.")
    check_file_exists(args.params_filename)

    if args.forcefield_path is not None:
        check_file_exists(args.forcefield_path+'/ffbonded.itp')
        check_file_exists(args.forcefield_path+'/ffnonbonded.itp')
    else:
        # print('The default 2016h66 force field will be used.')
        # root_path = os.path.dirname(__file__)
        # args.forcefield_path=root_path+'/gromos2016h66.ff'
        print("No force field was specified, default values for the interactions parameter will be used.")
        args.forcefield_path = None

    if args.stepLength <= 0:
        quitError("Please specify a positive value for the step lenght.")

    # The following can only happen, if the user specifically overrides the default with None,
    # which would be a very strange thing to do...:
    if args.ngen is None:
        quitError("Please specify the number of generations to be built by using --ngen option.")

    if args.nsteps is None:
        quitError("Please specify the numberof optimize geometry steps by using --nsteps option.")

    if args.ngenga is None:
        quitError("Please specify the number of GA generations by using --ngenga option.")

    if args.npop is None:
        quitError("Please specify the number oof GA population by using --npop option.")
    ##################


def process_buildingmode(args):
    building_mode = DENDRIMER
    if args.polymer and args.dendrimer:
        quitError("Choose EITHER dendrimer OR polymer")
    if args.polymer:
        building_mode = POLYMER
    sanity_check_buildingode(building_mode, args)
    return building_mode


def process_outputformat(args):
    output_format = GROMACS_FORMAT  # default
    if args.gromacs and args.gromos:
        quitError("Choose EITHER gromacs OR gromos")
    if args.gromos:
        output_format = GROMOS_FORMAT
    return output_format


def process_arguments(args):
    if args.debug:
        print("You used the DEBUG flag. Additional print messages will be made in the course of the processing...")
    ff_param = ff_params(args.params_filename)
    ff_param = give_unique_identifier_to_dihedral_params(ff_param)

    building_mode = process_buildingmode(args)
    output_format = process_outputformat(args)
    return ff_param, building_mode, output_format


def debug_print_args(args):
    print("Printing Input arguments...")
    for arg in vars(args):
        print("\t{}\t{}".format(arg, getattr(args, arg)))
    print("... arguments over.\n\n")


def disentangle_args(args):
    args = [args.output_filename, args.core_filename, args.ter_filename, args.inter_filename, args.gro_filename,
            args.ngenga, args.ngen, args.nsteps, args.npop,
            args.connections_paths, args.topName, args.bbs_paths,
            args.debug, args.params_filename, args.verbose, args.nogeom,
            args.forcefield_path, args.nskipLJ, args.stepLength]
    return args
