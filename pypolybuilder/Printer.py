'''
Created on Apr 20, 2015

@author: vitor
'''
from Utils import output_filename, gro_filename, topName
import Dihedral


class Printer(object):
    '''
    classdocs
    '''


    def __init__(self, params):
        '''
        Constructor
        '''

    @classmethod
    def printGromacs(cls, top):
        name = topName if topName else 'NEWTOP'
        with open(output_filename, 'w') as f:
            f.write("[ moleculetype ]\n")
            f.write("; Name   nrexcl\n")
            f.write(name + "    2\n\n")
            f.write("[ atoms ]\n")
            f.write(";  nr  type  resnr  resid  atom  cgnr  charge    mass    total_charge\n")

            for atom in top.get_atom_list():
                f.write("    " + str(atom.get_nr()) + "    " + str(atom.get_atomtype()) + "\t" + str(atom.get_resnr()) + "\t"
                    + str(atom.get_resid()) + "\t" + str(atom.get_atom()) + "\t" + str(atom.get_cgnr()) + "\t"
                    + str(atom.get_charge()) + "\t" + str(atom.get_mass()) + "\n")

            f.write("\n")
            f.write("[ bonds ]\n")
            f.write(";  ai   aj  funct   c0         c1 \n")

            for bond in top.get_bond_list():
                a1 = "%6s" % str(bond.get_a_1().get_nr())
                a2 = "%6s" % str(bond.get_a_2().get_nr())
                func = "%10s" % str(bond.get_func())
                param = "%10s" % str(bond.get_param())
                f.write(a1 + a2 + func + param + "\n")

            f.write("\n")
            f.write("[ angles ]\n")
            f.write(";  ai   aj   ak  funct   angle     fc\n")

            for angle in top.get_angle_list():
                a1 = "%6s" % str(angle.get_a_1().get_nr())
                a2 = "%6s" % str(angle.get_a_2().get_nr())
                a3 = "%6s" % str(angle.get_a_3().get_nr())
                func = "%10s" % str(angle.get_func())
                param = "%10s" % str(angle.get_param())
                f.write(a1 + a2 + a3 + func + param + "\n")

            f.write("\n")
            f.write("[ dihedrals ]\n")
            f.write(";  ai   aj   ak   al  funct    ph0      cp     mult\n")

            for dihedral in top.get_dihedral_list():
                a1 = "%6s" % str(dihedral.get_a_1().get_nr())
                a2 = "%6s" % str(dihedral.get_a_2().get_nr())
                a3 = "%6s" % str(dihedral.get_a_3().get_nr())
                a4 = "%6s" % str(dihedral.get_a_4().get_nr())
                func = "%10s" % str(dihedral.get_func())
                param = "%10s" % str(dihedral.get_param())
                f.write(a1 + a2 + a3 + a4 + func + param + "\n")

            f.write("\n")
            
            f.write("[ exclusions ]\n")
            for atom in top.get_atom_list():
                for j in atom.get_exclusion_extra().get_exclusion_extra():
                    f.write(str(atom.get_nr()) + "   " + str(j) + "\n")
                
        # f.close()

    @classmethod
    def printGROMACSCoordFile(cls,top,boxdim=[1.82060, 1.82060, 1.82060]):
        with open(gro_filename, 'w') as f:
            f.write("Positions\n")
            f.write("%5s" % str(len(top.get_atom_list())) + "\n")

            for atom in top.get_atom_list():
                a = "%5s" % atom.get_resnr()
                b = "%-5s" % atom.get_resid()
                c = "%5s" % atom.get_atom()
                d = "%5s" % atom.get_nr()
                x = "%8.3f" % atom.get_x()
                y = "%8.3f" % atom.get_y()
                z = "%8.3f" % atom.get_z()
                #x = 0
                #y = 0
                #z = 0
                f.write(str(a) + str(b) + str(c) + str(d) + str(x) + str(y) + str(z) + "\n")

            # f.write("   {:8}   {:8}   {:8}\n".format(*boxdim))
            f.write("   1.82060   1.82060   1.82060\n".format(*boxdim))
        # f.close()
    @classmethod
    def printGROMOSCoordFile(cls,top):
        with open(gro_filename, 'w') as f:
            f.write("TITLE\n        GromosXX\n        Automatically generated input file\n        hortab Sun May 12 17:49:41 2013")
            f.write("\n        final configuration\nEND\n")
            f.write("TIMESTEP\n        1000000     0.000000000\nEND\nPOSITION\n# first 24 chars ignored\n")
            #f.write(str(len(top.get_atom_list())) + "\n")

            for atom in top.get_atom_list():
                a = "%5s" % atom.get_resnr()
                b = "%-6s" % atom.get_resid()
                c = "%-5s" % atom.get_atom()
                d = "%7s" % atom.get_nr()
                x = "%15.9f" % atom.get_x()
                y = "%15.9f" % atom.get_y()
                z = "%15.9f" % atom.get_z()
                #x = 0
                #y = 0
                #z = 0
                f.write(str(a) + " " + str(b) + str(c) + str(d) + str(x) + str(y) + str(z) + "\n")

            f.write("END\nGENBOX\n    1\n    3.552498053    3.552498053    4.126560555\n")
            f.write("   90.000000000   90.000000000   90.000000000\n")
            f.write("    0.000000000    0.000000000    0.000000000\n")
            f.write("    0.000000000    0.000000000    0.000000000\n")
            f.write("END")
        # f.close()

    @classmethod
    def printGromos(cls,top):
        with open(output_filename, 'w') as f:
            f.write("MTBUILDBLSOLUTE\n#@BLOCKTYPE\n# building block (residue, nucleotide, etc.)\n# RNME\n")
            f.write("NEWTOP\n")
            f.write("# number of atoms, number of preceding exclusions\n# NMAT NLIN\n")
            f.write(str(len(top.get_atom_list())) + "    0\n# preceding exclusions\n#ATOM\n")
            f.write("# atoms\n#ATOM ANM  IACM MASS        CGMICGM MAE MSAE\n")

            for atom in top.get_atom_list():
                j1 = "%5s" % str(atom.get_nr())
                j2 = "%-5s" % str(atom.get_atom())
                j3 = "%4s" % str(atom.get_atomtype())
                j4 = "%10s" % str(atom.get_mass())
                j5 = "%11s" % str(atom.get_charge())
                j6 = "%4s" % str(atom.get_cgnr())
                j7 = "%4s" % str(len(atom.get_exclusion_list().get_exclusion_list()))
                j8 = ""
                for j in atom.get_exclusion_list().get_exclusion_list():
                    j8 += "%5s" % str(j)
                f.write(j1 + " " + j2 + j3 + j4 + j5 + j6  + j7 + j8 + "\n")

            f.write("# bonds\n#  NB\n")
            bl = "%5s" % str(len(top.get_bond_list()))
            f.write(bl + "\n")
            f.write("#  IB   JB  MCB\n")

            for bond in top.get_bond_list():
                a1 = "%5s" % str(bond.get_a_1().get_nr())
                a2 = "%5s" % str(bond.get_a_2().get_nr())
                #func = "%10s" % str(bond.get_func())
                param = "%5s" % str(bond.get_param()[3:])
                f.write(a1 + a2 + param + "\n")

            f.write("# bond angles\n#  NBA\n")
            al = "%5s" % str(len(top.get_angle_list()))
            f.write(al + "\n")
            f.write("#  IB   JB   KB  MCB\n")

            for angle in top.get_angle_list():
                a1 = "%5s" % str(angle.get_a_1().get_nr())
                a2 = "%5s" % str(angle.get_a_2().get_nr())
                a3 = "%5s" % str(angle.get_a_3().get_nr())
                #func = "%10s" % str(angle.get_func())
                param = "%5s" % str(angle.get_param()[3:])
                f.write(a1 + a2 + a3 + param + "\n")

            impropers = []
            dihedrals = []

            for dihedral in top.get_dihedral_list():
                if(dihedral.get_func() == 1):
                    dihedrals.append(dihedral)
                else:
                    impropers.append(dihedral)


            f.write("# improper dihedrals\n# NIDA\n")
            il = "%5s" % str(len(impropers))
            f.write(il + "\n")
            f.write("#  IB   JB   KB   LB  MCB\n")

            for dihedral in impropers:
                a1 = "%5s" % str(dihedral.get_a_1().get_nr())
                a2 = "%5s" % str(dihedral.get_a_2().get_nr())
                a3 = "%5s" % str(dihedral.get_a_3().get_nr())
                a4 = "%5s" % str(dihedral.get_a_4().get_nr())
                #func = "%10s" % str(dihedral.get_func())
                param = "%5s" % str(dihedral.get_param()[3:])
                f.write(a1 + a2 + a3 + a4 + param + "\n")


            f.write("# dihedrals\n# NDA\n")
            dl = "%5s" % str(len(dihedrals))
            f.write(dl + "\n")
            f.write("#  IB   JB   KB   LB  MCB\n")


            for dihedral in dihedrals:
                a1 = "%5s" % str(dihedral.get_a_1().get_nr())
                a2 = "%5s" % str(dihedral.get_a_2().get_nr())
                a3 = "%5s" % str(dihedral.get_a_3().get_nr())
                a4 = "%5s" % str(dihedral.get_a_4().get_nr())
                #func = "%10s" % str(dihedral.get_func())
                param = "%5s" % str(dihedral.get_param()[3:])
                f.write(a1 + a2 + a3 + a4 + param + "\n")

            f.write("# LJ exceptions\n# NEX\n    0\n#@FREELINE\nEND")
        # f.close()


    @classmethod
    def printTRAJFile(cls,top):
        with open("Traj-DEBUG.gro", 'w') as f:
            f.write("Positions\n")
            f.write("%5s" % str(len(top.get_atom_list())) + "\n")

            for atom in top.get_atom_list():
                a = "%5s" % atom.get_resnr()
                b = "%-5s" % atom.get_resid()
                c = "%5s" % atom.get_atom()
                d = "%5s" % atom.get_nr()
                x = "%8.3f" % atom.get_x()
                y = "%8.3f" % atom.get_y()
                z = "%8.3f" % atom.get_z()
                #x = 0
                #y = 0
                #z = 0
                f.write(str(a) + str(b) + str(c) + str(d) + str(x) + str(y) + str(z) + "\n")

            f.write("   1.82060   1.82060   1.82060\n")
        # f.close()
