'''
Created on Jan 23, 2015

@author: Vitor
'''
import argparse
from FFParams import ff_params
import getopt
import sys
from constants import DENDRIMER
from constants import POLYMER
from constants import GROMOS_FORMAT
from constants import GROMACS_FORMAT
import shutil
import sys
import os.path
import numpy as np

script_dir = os.path.dirname(__file__)

logoFilename = script_dir + '/asciilogo.txt'
with open(logoFilename, "r") as f:
    shutil.copyfileobj(f, sys.stdout)

output_filename  = 'default.itp'
core_filename    = None
ter_filename     = None
inter_filename   = None
gro_filename     = 'default.gro'
params_filename  = None
topName          = None
connections_paths = None
optimizeGeometry = False
bbs_paths = []
buildingMode = DENDRIMER
outputFormat = GROMACS_FORMAT
ngen     = 0 #number of generations

parser = argparse.ArgumentParser(description='All available options are shown below')
parser.add_argument('--core', help='Core filename')
parser.add_argument('--inter', help='Inter filename')
parser.add_argument('--ter', help='Ter filename')
parser.add_argument('--output', help='Output filename')
parser.add_argument('--params', help='Forcefield params filename')
parser.add_argument('--ngen', help='Number of generations to be built')
parser.add_argument('--optgeo',action='store_true',  help='If enabled pypolybuilder will try to find a reasonable geometry')
parser.add_argument('--polymer',action='store_true', help='If enabled pypolybuilder will expect Polymer inputs and build it')
parser.add_argument('--bbs', help='Building blocks filename used to build a Polymer')
parser.add_argument('--in', help='Building blocks connections filename')
parser.add_argument('--gromos', action='store_true',help='If enabled pypolybuilder will print output with gromos format')
parser.add_argument('--name', help='Topology name to be printed on output file')

args = parser.parse_args()

options, remainder = getopt.getopt(sys.argv[1:], 'o:c', ['output=',
                                                         'core=',
                                                         'inter=',
                                                         'ter=',
                                                         'gro=',
                                                         'params=',
                                                         'ngen=',
                                                         'in=',
                                                         'optgeo',
                                                         'polymer',
                                                         'gromos',
                                                         'name=',
                                                         'bbs=',
                                                         'h'
                                                         ])


for opt, arg in options:
    if opt in ('-o', '--output'):
        output_filename = arg
    elif opt in ('-c', '--core'):
        core_filename = arg
    elif opt in ('-t', '--ter'):
        ter_filename = arg
    elif opt in ('-i', '--inter'):
        inter_filename = arg
    elif opt in ('-g', '--gro'):
        print("found gfname")
        gro_filename = arg
    elif opt in ('-l', '--params'):
        params_filename = arg
    elif opt in ('-n','--ngen'):
        ngen = int(arg);
    elif opt in ('-in','--in'):
        connections_paths = arg
    elif opt in ('-optgeo', '--optgeo'):
        optimizeGeometry = True
    elif opt in ('-polymer', '--polymer'):
        buildingMode = POLYMER
    elif opt =='-gromos' or opt=='--gromos':
        print("found gromos")
        outputFormat = GROMOS_FORMAT
    elif opt in ('-name', '--name'):
        topName = arg
    elif opt in ('-bbs','--bbs'):
        arg = arg.split(',')
        for bb_filename in arg:
            bbs_paths.append(bb_filename)

ff_param = ff_params(params_filename)

def quitError(message):
    print("An error occurred.")
    print(message)
    print("See pypolybuilder -h for more instructions.")
    exit()

if(params_filename == None or os.path.exists(params_filename) == False):
    quitError("Please specify an existing Force Field parameters filename by using --params option.")

if(ngen == None):
    quitError("Please specify the number of generations to be built by using --ngen option.")

if(buildingMode == POLYMER):
    if(len(bbs_paths) < 1):
        quitError("Please specify an existing  polymer building blocks by using --bbs option.")
    if (connections_paths == None or os.path.exists(connections_paths) == False):
        quitError("Please specify an existing  polymer connections filename by using --in option.")

if(buildingMode == DENDRIMER):
    if(core_filename == None or os.path.exists(core_filename) == False):
        quitError("Please specify an existing  dendrimer core filename by using --core option.")
    if (ter_filename == None or os.path.exists(ter_filename) == False):
        quitError("Please specify an existing  dendrimer terminal filename by using --ter option.")


def unit_vector(vector):
    return vector / np.linalg.norm(vector)


def vectorFromAtoms(atom1, atom2):
    vec12 = [atom1.get_x() - atom2.get_x(), atom1.get_y() - atom2.get_y(),
             atom1.get_z() - atom2.get_z()]
    return vec12