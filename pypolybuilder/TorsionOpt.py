'''
Created on Jan 22, 2019

@author: bruno
Edit...
@author: paquoika
'''

from __future__ import division


from Printer import Printer
# import sys
import time
import numpy as np
import math
import copy
import Dihedral
from Polymer import Polymer
from Utils import unit_vector, deg_to_rad
from Utils import vectorFromAtoms
from vbga import VBGA, calcSphericity, RoulettSelectionMethod, AritmeticCrossoverMethod, DefaultMutation
from scipy.spatial import distance
from Utils import ngenga, npop
from Utils import debug as DEBUG, print_pypoly_warning
from Utils import verbose as args_verbose
from pypoly_io import quitError


Matrix = None
# HARD-CODED PARAMETERS FOR GENERATING INITIAL STRUCTURE
bonddist = [1, 1, 1.09, 1.12, 1.23, 1.25, 1.32, 1.33, 1.33, 1.33, 1.34, 1.34, 1.36, 1.38, 1.39, 1.39, 1.4, 1.43,
            1.43, 1.435, 1.47, 1.48, 1.48, 1.48, 1.5, 1.52, 1.53, 1.61, 1.63, 1.78, 1.78, 1.83, 1.87, 1.98, 2,
            2.04, 2.21, 1, 1.1, 1.758, 1.53, 1.93799, 1.76, 1.265, 1.35, 1.63299, 2.33839, 2.90283, 2.79388,
            2.91189, 2.077, 2.87407]
bondangle = [90.000000, 90.000000, 96.000000, 100.000000, 103.000000, 104.000000, 108.000000, 109.500000,
             109.500000, 109.500000, 109.500000, 109.500000, 109.500000, 109.600000, 111.000000, 113.000000,
             115.000000, 115.000000, 115.000000, 116.000000, 116.000000, 117.000000, 120.000000, 120.000000,
             120.000000, 120.000000, 120.000000, 120.000000, 120.000000, 121.000000, 122.000000, 123.000000,
             124.000000, 125.000000, 125.000000, 126.000000, 126.000000, 126.000000, 132.000000, 155.000000,
             180.000000, 109.500000, 107.570000, 111.300000, 97.400000, 106.750000, 108.530000, 109.500000,
             107.600000, 109.500000, 110.300000, 111.400000, 117.200000, 121.400000]
improper = [1000, 0.0, 35.26439, 0.0, 180.0, -35.26439]


# STATIC FUNCTIONS --> may be put to utils
def sanity_check_identical_vals(a, k, ZangleTo_out):
    if k != a or ZangleTo_out != a:
        raise ValueError("{} {} {}".format(a, k, ZangleTo_out))


def get_respective_bond_angle_param(angle):
    index_angle_param = int(angle.get_param()[3:])
    return bondangle[index_angle_param]  # bond angle is a global variable list


def print_dihedral(dihed):
    alist = [dihed.get_a_1(), dihed.get_a_2(), dihed.get_a_3(), dihed.get_a_4()]
    rep = ["{}{}".format(atom.get_atom(), atom.get_nr()) for atom in alist]
    print("\t".join(rep))


def normalize_vector(v, catch0=True):
    vnorm = np.linalg.norm(v)
    if catch0 and vnorm < 1e-5:
        vnorm = 1e-5
    return v/vnorm


def flag_members_of_tetraedrals(atom):
    alist = atom.get_atom_list()
    for atom in alist:
        atom.set_is_neighbor_of_tetraedral(True)


class Zmatrix(object):
    atoms = []
    boundto = []
    bondvalue = []
    angleto = []
    angvalue = []
    dihedralto = []
    dihvalue = []

    def __init__(self, atoms, boundto, bond, angleto, ang, dihedralto, dih):
        self.atoms = atoms
        self.boundto = boundto
        self.bondvalue = bond
        self.angleto = angleto
        self.angvalue = ang
        self.dihedralto = dihedralto
        self.dihvalue = dih

    def sanity_check(self):
        # No atomID should negative
        temp = np.array([atom.get_nr() for atom in self.atoms]).flatten()
        if True in (temp<0):
            self.printDebugMatrix()
            raise ValueError("Something is wrong with your Zmatrix:\nNegative atomID")

        temp = np.array([bond for bond in self.boundto]).flatten()
        if True in (temp<0):
            self.printDebugMatrix()
            raise ValueError("Something is wrong with your Zmatrix:\nNegative atomID of bound partner")

        temp = np.array([bond for bond in self.angleto]).flatten()
        if True in (temp < 0):
            self.printDebugMatrix()
            raise ValueError("Something is wrong with your Zmatrix:\nNegative atomID of angle partner")

        temp = np.array([bond for bond in self.dihedralto]).flatten()
        if True in (temp < 0):
            self.printDebugMatrix()
            raise ValueError("Something is wrong with your Zmatrix:\nNegative atomID of dihedral partner")
        # bonds should be finite
        temp = np.array([bond for bond in self.bondvalue]).flatten()
        if (True in (temp < 0)) or (False in np.isfinite(temp)):
            self.printDebugMatrix()
            raise ValueError("Something is wrong with your Zmatrix:\nNonfinite or negative bond length")
        # angles should be finite
        temp = np.array([angle for angle in self.angvalue]).flatten()
        if False in np.isfinite(temp):
            self.printDebugMatrix()
            raise ValueError("Something is wrong with your Zmatrix:\nNonfinite angle")
        return True

    def printMatrix(self):
        natoms = len(self.atoms)
        print(natoms)
        print('{:>4s}'.format(self.atoms[0].get_atomtype()[0]))
        print('{:>4s} {:>4d} {:>11.5f}'.format(self.atoms[1].get_atomtype()[0],
                                               self.boundto[0], float(self.bondvalue[0])))
        print('{:>4s} {:>4d} {:>11.5f} {:>4d} {:>11.5f}'.format(self.atoms[2].get_atomtype()[0],
                                                                self.boundto[1], float(self.bondvalue[1]),
                                                                self.angleto[0], float(self.angvalue[0])))
        for i in range(3, natoms):
            print('{:>4s} {:>4d} {:>11.5f} {:>4d} '
                  '{:>11.5f} {:>4d} {:>11s}'.format(self.atoms[i].get_atom(), self.boundto[i - 1],
                                                    float(self.bondvalue[i - 1]), self.angleto[i - 2],
                                                    float(self.angvalue[i - 2]), self.dihedralto[i - 3],
                                                    str(self.dihvalue[i - 3])))

    def printDebugMatrix(self):
        natoms = len(self.atoms)
        print("Number of atoms: {}".format(natoms))
        # table header:
        print('{:>4} {:>4} {:>8} {:>4} {:>8} {:>8} {:>6} {:>6}'.format("id", "type", "boundto", "b-dist",
                                                                       "angleTo", "ang", "DihTo", "Dih"))
        # first three lines are special, bc there are no bonds/ angles/ dihedrals
        print('{:>4d} {:>4s}'.format(1, self.atoms[0].get_atomtype()[0]))
        print('{:>4d} {:>4s} {:>4d} {:>11.5f}'.format(2, self.atoms[1].get_atomtype()[0],
                                                      self.boundto[0], float(self.bondvalue[0])))
        print('{:>4d} {:>4s} {:>4d} {:>11.5f} '
              '{:>4d} {:>11.5f}'.format(3, self.atoms[2].get_atomtype()[0], self.boundto[1],
                                        float(self.bondvalue[1]), self.angleto[0], float(self.angvalue[0])))
        for i in range(3, natoms):
            try:
                try: bt = self.boundto[i - 1]
                except IndexError:
                    bt = -1
                try: bv = float(self.bondvalue[i - 1])
                except IndexError:
                    bv = np.nan
                try: at = self.angleto[i - 2]
                except IndexError:
                    at = -1
                try: av = float(self.angvalue[i - 2])
                except IndexError:
                    av = np.nan
                try: dt = self.dihedralto[i - 3]
                except IndexError:
                    dt = -1
                try: dv = str(self.dihvalue[i - 3])
                except IndexError:
                    dv = str(np.nan)
                print('{:>4d} {:>4s} {:>4d} '
                      '{:>11.5f} {:>4d} {:>11.5f} {:>4d} {:>11s}'.format(i+1, self.atoms[i].get_atom(),
                                                                         bt, bv, at, av, dt, dv))
            except IndexError:
                print("IndexError: failed at atom-id: {}".format(i+1))

    def printDihedralTo(self):
        for i in self.dihedralto:
            print(i)

    def randomize(self):
        for i in range(0, len(self.dihvalue)):
            if self.dihvalue[i] == "NULL":
                a = np.random.randint(12)
                newangle = float(a * 30)
                self.dihvalue[i] = newangle

    def movableTorsions(self):
        # Improve this at a later stage to include restrictions (maybe need info about dihedrals)
        indexList = []
        for i in range(0, len(self.dihvalue)):
            if self.dihvalue[i] == "NULL":
                indexList.append(i)
        return indexList

    def converToXYZ(self):
        natoms = len(self.atoms)

        xyzarr = np.zeros([natoms, 3])
        if natoms > 1:
            xyzarr[1] = [self.bondvalue[0], 0.0, 0.0]

        if natoms > 2:
            i = self.boundto[1] - 1
            j = self.angleto[0] - 1
            r = self.bondvalue[1]
            theta = deg_to_rad(self.angvalue[0])
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            a_i = xyzarr[i]
            b_ij = xyzarr[j] - xyzarr[i]
            if b_ij[0] < 0:
                x = a_i[0] - x
                y = a_i[1] - y
            else:
                x = a_i[0] + x
                y = a_i[1] + y
            xyzarr[2] = [x, y, 0.0]

        for n in range(3, natoms):
            # TODO put coordinate transform in separate function (and separate module)
            # Transform spherical coordinates to cartesian)
            r = self.bondvalue[n - 1]
            theta = deg_to_rad(self.angvalue[n - 2])
            phi = deg_to_rad(self.dihvalue[n - 3])

            sinTheta = np.sin(theta)
            cosTheta = np.cos(theta)
            sinPhi = np.sin(phi)
            cosPhi = np.cos(phi)

            x = r * cosTheta
            y = r * cosPhi * sinTheta
            z = r * sinPhi * sinTheta

            i = self.boundto[n - 1] - 1
            j = self.angleto[n - 2] - 1
            k = self.dihedralto[n - 3] - 1
            a = xyzarr[k]
            b = xyzarr[j]
            c = xyzarr[i]

            ab = b - a
            bc = c - b
            bc = normalize_vector(bc)
            nv = np.cross(ab, bc)
            nv = normalize_vector(nv)
            ncbc = np.cross(nv, bc)

            new_x = c[0] - bc[0] * x + ncbc[0] * y + nv[0] * z
            new_y = c[1] - bc[1] * x + ncbc[1] * y + nv[1] * z
            new_z = c[2] - bc[2] * x + ncbc[2] * y + nv[2] * z
            xyzarr[n] = [new_x, new_y, new_z]

        return (xyzarr)


class TorsionOpt(object):
    def __init__(self, top):
        self.top = top

    def topToZmatrix(self):
        # TODO Is this still used somewhere?
        global debug
        debug = False

        # copy lists, which will be popped below
        # print("copy bonds")  # debug
        list_bonds = copy.deepcopy(self.top.get_bond_list())
        # print("copy angles")  # debug
        list_angles = copy.deepcopy(self.top.get_angle_list())
        # print("copy dihedrals")  # debug
        list_dihedrals = copy.deepcopy(self.top.get_dihedral_list())
        if debug:
            print("DEBUG", np.shape(list_angles))

        # Check for tetrahedral centers
        for atom in self.top.get_atom_list():
            if self.top.is_tetracoordinated(atom):
                if args_verbose:
                    print("Atom " + str(atom.get_nr()) + " is a tetrahedral center")
                flag_members_of_tetraedrals(atom)
                list_dihedrals = self.createAllImpropersSorroundingThCenter(list_angles, list_dihedrals, atom)

        ##################################
        # Create Z-matrix
        Zatom = []
        ZbondTo, ZbondValue = [], []
        ZangleTo, ZangleValue = [], []
        ZdihedralTo, ZdihValue = [], []
        atom_list = self.top.atom_list
        for atom in atom_list:
            atom_id0 = int(atom.get_nr())
            Zatom.append(atom)  # Add atom
            # special case: first entry in the Z-matrix:
            if atom_id0 == 1:
                continue  # doesnt have bonds, angles, or dihedrals
            for index_bond, bond in enumerate(list_bonds):
                row_complied = False
                atom_id1 = int(bond.get_a_1().get_nr())  # Index bond of atom 1
                if atom_id0 == int(bond.get_a_2().get_nr()):  # Checks if the bond matches this atom
                    if debug:
                        print('DEBUG: BOND {} {} '
                              '[{}, {}]'.format(atom_id0,atom_id1,bond.get_a_1().get_nr(),bond.get_a_2().get_nr()))

                    ZbondTo.append(atom_id1)  # Add bond
                    index_bond_param = int(bond.get_param()[3:])
                    ZbondValue.append(bonddist[index_bond_param] / 10)  # Add bond value
                    list_bonds.pop(index_bond)  # Delete bond from lookup list
                    if atom_id0 == 3:  # Special case, there is no dihedral, but an angle here in the Z-matrix
                        # list_angles is being updated
                        temp_ZangleTo, temp_ZangleValue, list_angles = self.caseWithoutDihedral(list_angles,
                                                                                                atom_id0, atom_id1)
                        ZangleTo.append(temp_ZangleTo)
                        ZangleValue.append(temp_ZangleValue)
                    elif atom_id0 > 3:  # For all cases below, there is a dihedral
                        temp_return, updated_lists = self.caseWithDihedral(list_angles, atom_id0, atom_id1,
                                                                           row_complied, list_dihedrals, atom)
                        # update lists
                        list_angles, list_dihedrals = updated_lists
                        # append to Z-matrix
                        temp_ZangleTo, temp_ZangleValue, temp_ZdihedralTo, temp_ZdihValue = temp_return
                        ZdihedralTo.append(temp_ZdihedralTo)
                        ZdihValue.append(temp_ZdihValue)
                        ZangleTo.append(temp_ZangleTo)
                        ZangleValue.append(temp_ZangleValue)
                    # else: atomid 0, or 1, there is no angles or dihedrals for these atoms

        return Zatom, ZbondTo, ZbondValue, ZangleTo, ZangleValue, ZdihedralTo, ZdihValue

    def process_angle(self, atoms, atomid0, atomid1, list_angles, angle):
        """This routine is the same for Case with and Case Without dihedral."""
        ZangleTo_out, ZangleValue_out = -2, np.nan  # dummy values, if these appear, your Z-matrix will be pathological
        k=-1
        a1, a2, a3 = atoms
        if a3 == atomid0:
            if a2 == atomid1: #MAYK and a2 > a1:  # Type angle (9 - 8 - 7)
                # print(" a3 -> a2")
                k = int(angle.get_a_1().get_nr())
                ZangleTo_out = a1  # Add angle
                sanity_check_identical_vals(a1, k, ZangleTo_out)
                ZangleValue_out = get_respective_bond_angle_param(angle)

            elif a1 == atomid1:  #MAYK and a2 < a1:  # Type angle (9 - 7 - 8)
                # print(" a3 -> a1")
                k = int(angle.get_a_2().get_nr())
                ZangleTo_out = a2  # Add angle
                sanity_check_identical_vals(a2, k, ZangleTo_out)
                ZangleValue_out = get_respective_bond_angle_param(angle)

        elif a1 == atomid0:
            if a2 == atomid1:
                # print(" a1 -> a2")
                k = int(angle.get_a_3().get_nr())
                ZangleTo_out = a3  # Add angle
                sanity_check_identical_vals(a3, k, ZangleTo_out)
                ZangleValue_out = get_respective_bond_angle_param(angle)
            # Probably it will be needed in some future case
            # elif a3 == atomid1:
            #     # print(" a1 -> a2")
            #     k = int(angle.get_a_2().get_nr())
            #     ZangleTo_out = a2  # Add angle
            #     sanity_check_identical_vals(a2, k, ZangleTo_out)
            #     ZangleValue_out = get_respective_bond_angle_param(angle)

        elif a2 == atomid0:
            if a1 == atomid1:
                # print(" a2 -> a1")
                k = int(angle.get_a_3().get_nr())
                ZangleTo_out = a3  # Add angle
                sanity_check_identical_vals(a3, k, ZangleTo_out)
                ZangleValue_out = get_respective_bond_angle_param(angle)
            # elif a3 == atomid1:
            #     # print(" a2 -> a3")
            #     k = int(angle.get_a_1().get_nr())
            #     ZangleTo_out = a1  # Add angle
            #     sanity_check_identical_vals(a1, k, ZangleTo_out)
            #     ZangleValue_out = get_respective_bond_angle_param(angle)

        else:
            print("AtomID0 is not in the angle_tuple. This should never happen!!?!!")
        return ZangleTo_out, ZangleValue_out, list_angles, k

    def caseWithoutDihedral(self, list_angles, atomid0, atomid1):
        """Note: In here we do not pop the angles_list. The angle
        we find here, will stay in the lookup list.
        PAQ: I dont see a problem with that...
        """
        # TODO make zero and 1 counting consistent
        ZangleTo_out, ZangleValue_out = -3, np.nan  # dummy vals for debug
        for index_angle, angle in enumerate(list_angles):
            atoms = angle.get_all_atoms_nr()
            if atomid0 in atoms and atomid1 in atoms:
                ZangleTo_out, ZangleValue_out, list_angles, k = self.process_angle(atoms,atomid0,atomid1,
                                                                                   list_angles, angle)
                break  # found an abitrary angle, stop searching.

            else:  # didnt find any angles with these atoms. check next angle in list
                continue
        if set(np.sort([atomid0, atomid1, k])) != set(np.sort(atoms)):
            print('WARNING: ANGLE woDih NEQ {} {} {}; {} {} {}'.format(atomid0, atomid1, k,*atoms))
            vals = np.concatenate([np.sort([atomid0, atomid1, k]), np.sort(atoms)])
            print('WARNING: \t\t\t {} {} {}; {} {} {}'.format(*vals))
        return ZangleTo_out, ZangleValue_out, list_angles

    def get_dihed_ref_indices(self, dihedral):
        v = dihedral.get_all_atoms_nr()  # Vector containing the number (nr) of the atoms in a dihedral
        len_v = len(v)

        p_i = v.index(max(v))  # Index of the atom with largest number (nr)
        i_ref = v[p_i]  # this is the atom with max index in this vector
        if dihedral.is_torsional():  # if this is a torsional dihedral
            value = 'NULL'
            if p_i == 3:  # Normal dihedral with reverse grow
                j_ref = v[(p_i - 1)]
                k_ref = v[(p_i - 2)]
                l_ref = v[(p_i - 3)]
            else:  # Normal dihedral with direct grow
                j_ref = v[(p_i + 1) % len_v]
                k_ref = v[(p_i + 2) % len_v]
                l_ref = v[(p_i + 3) % len_v]
                # Shall we also consider the cases where p_i = 1 or p_i = 2?

        else:  # improper dihedral
            param = dihedral.get_param()
            try:
                index_dihedral_param = int(param.split("_")[1])  # split gd_***
                splitted_ref_idx = True
            except:  # This exception is extremely broad. So ultimately if there is ANYTHING WRONG
                # with the parameter this will catch
                warn_msg = "\nCouldn't split and/or convert dihedral parameter in your input >{}< " \
                           "according to GROMOS convention!!!\n" \
                           "For the optimization of the torsional degrees of freedom, " \
                           "PyPolyBuilder will try to make this improper dihedral a torsional dihedral. " \
                           "(Which is WILDLY unphysical, " \
                           "but a pragmatic solution to still get some initial structure!)\n" \
                           "This may mess up the optimized structure! " \
                           "Also make sure, that your topology is correct...\n".format(param)
                print_pypoly_warning(warn_msg, category=RuntimeWarning)
                splitted_ref_idx = False
            if splitted_ref_idx:
                try:
                    value = improper[index_dihedral_param]
                except IndexError:
                    print("IndexError after splitting {}, tried {} in list:".format(param, index_dihedral_param),
                          improper)
                    quitError("Problem in get_dihed_ref_indices...\n"
                              "Trying to access non-existent value in reference table "
                              "for torsions of improper dihedrals.")
            else:
                value = "NULL"

            j_ref = v[(p_i + 1) % len_v]
            k_ref = v[(p_i + 2) % len_v]
            l_ref = v[(p_i + 3) % len_v]

        return i_ref, j_ref, k_ref, l_ref, value

    def caseWithDihedral(self, list_angles, atomid0, atomid1, row_complied, list_dihedrals, atom, rollback=False):
        ZangleTo_out, ZangleValue_out = -3, np.nan  # dummy vals for debug
        ZdihedralTo_out, ZdihValue_out = -4, np.nan  # dummy vals for debug
        k = -1  # this default value WOULD be caught in the sanity_check
        for index_angle, angle in enumerate(list_angles):
            if row_complied:
                break  # Ends if line is complete
            atoms = angle.get_all_atoms_nr()  # Angle: triplet of atom numbers
            if atomid0 in atoms and atomid1 in atoms:
                ZangleTo_out, ZangleValue_out, list_angles, k = self.process_angle(atoms, atomid0, atomid1,
                                                                                   list_angles, angle)
            else:  # didnt find any angles with these atoms. check next angle in list
                continue
            if set(np.sort([atomid0, atomid1, k])) != set(np.sort(atoms)):
                print('WARNING: This should never happen: Dihedral, should have had these indices:\n'
                      ' {} {} {}; {} {} {}'.format(atomid0, atomid1, k, *atoms))
                vals = np.concatenate([np.sort([atomid0, atomid1, k]), np.sort(atoms)])
                print('WARNING: \t\t\t {} {} {}; {} {} {}'.format(*vals))
            dihedral_possible = []  # Make list with possible dihedrals
            dihedral_possible_value = []
            for index_dihedral, dihedral in enumerate(list_dihedrals):
                i_ref, j_ref, k_ref, l_ref, value = self.get_dihed_ref_indices(dihedral)

                if atom.get_is_neighbor_of_tetraedral():  # atom is in tetrahedral center
                    pass
                # else:  # regular operation:
                if i_ref == atomid0 and j_ref == atomid1 and k_ref == k:  # Dihedral compatibility test
                    l = l_ref
                    dihedral_possible.append(l)  # Save all possible dihedrals
                    dihedral_possible_value.append(value)
                    row_complied = True  # Signal for complete row

                    if not atom.get_is_neighbor_of_tetraedral():  # atom is not around tetrahedral center
                        list_dihedrals.pop(index_dihedral)  # Dihedrals are never needed twice. is that correct?
                        # except for at tetrahedral centers...
            if len(dihedral_possible) > 0:
                ZdihedralTo_out = min(dihedral_possible)  # Choose the lower index
                ZdihValue_out = dihedral_possible_value[dihedral_possible.index(min(dihedral_possible))]
                # pop the angle list only if dihedral found
                if not atom.get_is_neighbor_of_tetraedral():  # atom is not around tetrahedral center
                    list_angles.pop(index_angle)  # never need angles, twice, if I found the dihedral behind it...?
                    # except for at tetrahedral centers...

            # else:  # no possible dihed "
        if ZdihedralTo_out == -4:
            print("############################## DIDNT FIND DIHEDRAL:")
            print(atomid0, atomid1, k)
            for dihedral in list_dihedrals:
                i_ref, j_ref, k_ref, l_ref, value = self.get_dihed_ref_indices(dihedral)
                if i_ref == atomid0 and j_ref == atomid1:
                    print_dihedral(dihedral)
            # if not rollback:
            #     print("Rollback and try with different angle...")
            #     list_angles_rollback = copy.deepcopy(list_angles)
            #     print("Removing angle from lookup table: {} : {}-{}-{}".format(index_angle, atomid0, atomid1, k))
            #     list_angles_rollback.pop(index_angle)
            #     # This function calls itself with a reduced list of angles
            #     roll_back_out = self.caseWithDihedral(list_angles_rollback, atomid0, atomid1,
            #                                           row_complied, list_dihedrals, atom, rollback=True)
            #     (ZangleTo_out, ZangleValue_out, ZdihedralTo_out, ZdihValue_out), \
            #     (list_angles_rollback, list_dihedrals) = roll_back_out

        return (ZangleTo_out, ZangleValue_out, ZdihedralTo_out, ZdihValue_out), (list_angles, list_dihedrals)

    def createAllImpropersSorroundingThCenter(self, list_angles, list_dihedrals, aT):
        # TODO list_angles is actually not used here
        atoms_of_Th_center = self.top.get_atoms_bound_to_me(aT)
        new_dihedrals = []
        for i in range(0, 4):
            for j in range(0, 4):
                # if self.top.triplet_belongs_to_angle(atoms_of_Th_center[i], atoms_of_Th_center[j], aT):
                if atoms_of_Th_center[i].get_nr() > atoms_of_Th_center[j].get_nr():
                    atomi = atoms_of_Th_center[i]
                    atomj = atoms_of_Th_center[j]
                    atomk = max(atoms_of_Th_center, key=lambda x: x.get_nr() if x.get_nr() < atomi.get_nr() else 0)
                    if atomk != atomj:
                        new_dihedral = Dihedral.Dihedral(atomi, aT, atomj, atomk, 2, "gi_2")
                        new_dihedrals.append(new_dihedral)
                        break
        if DEBUG:
            print("Created the following improper dihedrals at improper dihedral center :")
        for new_dihedral in new_dihedrals:
            if DEBUG:
                print(new_dihedral)
            list_dihedrals.append(new_dihedral)
        if args_verbose:
            print("\n")
        return list_dihedrals

    def addXYZtoTOP(self, xyz):
        n = 0
        for atom in self.top.get_atom_list():
            atom.set_x(xyz[n][0])
            atom.set_y(xyz[n][1])
            atom.set_z(xyz[n][2])
            n = n + 1

    def run(self):
        start = time.time()
        print("Converting Top to Z-matrix ...\n")
        a, b, c, d, e, f, g = self.topToZmatrix()
        end = time.time()
        print("Time for Z-Matrix construction: {:.3f} seconds.".format(end - start))
        global Matrix
        Matrix = Zmatrix(a, b, c, d, e, f, g)
        if Matrix.sanity_check():  # will raise en error, if any sanity checks is not passed
            print("Z-matrix successfully passed basic sanity checks.")
        if DEBUG:
            Matrix.printDebugMatrix()

        myGA = VBGA(Individual, calcSphericity, RoulettSelectionMethod, AritmeticCrossoverMethod, crossoverRate=80.0,
                    mutationMethod=DefaultMutation, mutationRate=20.0, populationSize=npop)

        print("\nPerforming optimization by genetic algorithm: "
              "Generation number = {}; Population = {}\nThis may take a while...\n".format(ngenga,npop))
        myGA.run(ngenga)

        xyz = myGA.getBest().matrix.converToXYZ()
        self.addXYZtoTOP(xyz)


class Individual:

    def __init__(self):
        self.matrix = copy.deepcopy(Matrix)
        self.imp = []
        self.fitValue = 0

    def randomize(self):
        for i in range(0, len(self.matrix.dihvalue)):
            if self.matrix.dihvalue[i] == "NULL":
                a = np.random.randint(24)
                newangle = float(a * 15)
                self.matrix.dihvalue[i] = newangle
                self.imp.append(1)
            else:
                self.imp.append(0)
