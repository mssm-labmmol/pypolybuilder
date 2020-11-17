'''
Created on Dec 22, 2013

@author: root
'''

from Itp      import Itp
from Angle    import Angle
from Bond     import Bond
from Dihedral import Dihedral
import random
import copy
import Utils
from Utils import ngen
from Utils import connections_paths
from Utils import bbs_paths
from BuildingBlockFactory import BuildingBlockFactory
from Connection import Connection
import numpy as np


def is_comment(line, comment_chars=[';', '#']):
    return line[0] in comment_chars


def lines_have_keyword(lines, keyword="[ BUILDING BLOCKS ]"):
    return np.sum([(keyword in line) for line in lines])


class Polymer(Itp):

    def create_coords(self):
        for atom in self.get_atom_list():
            x = random.uniform(0, 1)
            y = random.uniform(0, 1)
            z = random.uniform(0, 1)
            atom.set_x(x)
            atom.set_y(y)
            atom.set_z(z)

    def read_building_blocks(self, lines, lidx=0):
        bb_lines = []
        for line in lines[lidx:]:
            if is_comment(line):
                continue
            if line.isspace():
                break  # first blank line will stop reading after keyword
            if line[0] != ";" and "[ CONNECTS ]" not in line:
                bb_lines.append(line.split())
            else:
                break
        return bb_lines

    def process_building_blocks(self,bb_lines):
        for line in bb_lines:
            bb_number = int(line[0])
            bb_resname = line[1]
            building_block = self.bb_factories[bb_resname].makeBBForItp(self)
            self.building_blocks[bb_number] = building_block
            self.connect_bb(building_block)
            self.dict_number_resname[bb_number] = bb_resname

    def read_connects(self, lines, lidx=0):
        con_lines = []
        for line in lines[lidx:]:  # no connections should be given before the connection keyword
            if is_comment(line):
                continue
            if line.isspace():
                break  # first blank line will stop reading after keyword
            if not is_comment(line):
                line = line.split()
                if len(line) > 0:
                    if len(line) == 4:
                        con_lines.append(line)
                    elif len(line) < 4:
                        print("Error in connection file.")
                        quit()

        return con_lines

    def read_connects_w_dihedral(self, lines, lidx=0):
        con_lines = []
        for line in lines[lidx:]:  # no connections should be given before the connection keyword
            if is_comment(line):
                continue
            if line.isspace():
                break  # first blank line will stop reading after keyword
            if not is_comment(line):
                line = line.split()
                if len(line) > 0:
                    # print(line)
                    if len(line) == 8:
                        con_lines.append(line)
                    else:
                        print("Error in connection file. "
                              "Suspicious input for connections with dihedral:\n{}".format(line))
                        quit()

        return con_lines

    def preprocess_diheds_at_connects(self, con_lines_w_dihedral):
        if len(con_lines_w_dihedral) > 0:
            block_idxs = np.array([line[:2] for line in con_lines_w_dihedral], dtype=int)
            atom_idxs = np.array([line[2:6] for line in con_lines_w_dihedral], dtype=int)

            dihedral_specs = np.array([line[2:] for line in con_lines_w_dihedral], dtype=str)
            connect_ref_vals = np.array(list(zip(block_idxs[:, 0], block_idxs[:, 1], atom_idxs[:, 1], atom_idxs[:, 2])))
            return {"connect_ref_vals": connect_ref_vals,
                    "dihedral_specs": dihedral_specs}
        else:
            return {"connect_ref_vals": np.full((4, 1), np.nan),
                                     "dihedral_specs": np.full((4, 1), -1)}

    def get_specified_dihedral_evtl(self, con_line, specified_connections):
        # BBI    BBJ   IAT-1   IAT   JAT   JAT+1          funct    ph0      cp     mult
        # defaults, if none is found
        connect_vals = con_line[0], con_line[1], con_line[2], con_line[3]
        line = np.array(con_line, dtype=int)
        has_specified_dihedral = False
        dihedral_spec = []

        for vals, specs in zip(specified_connections["connect_ref_vals"],
                               specified_connections["dihedral_specs"]):
            # print(line, vals)
            # print(line[0],vals[0])
            if (line[0] == vals[0]) and (line[1] == vals[1]) and \
                    (line[2] == vals[2]) and (line[3] == vals[3]):
                has_specified_dihedral = True
                # print("FOUND")
                dihedral_spec.append(specs)
        dihedral_spec = {"atom_nrs": [int(x) for x in specified_connections["dihedral_specs"][0][:4]],
                         "funct": [int(x[4]) for x in dihedral_spec],
                         "ph0": [x[5] for x in dihedral_spec]}
        return has_specified_dihedral, connect_vals, dihedral_spec

    def process_connects(self, con_lines, con_lines_w_dihedral=[]):
        # print(con_lines_w_dihedral)

        specified_connections = self.preprocess_diheds_at_connects(con_lines_w_dihedral)

        for con_line in con_lines:  # iterate over connections
            has_spec_dihedral, connect_vals, dihedral_spec = self.get_specified_dihedral_evtl(con_line,
                                                                                              specified_connections)
            # print(connect_vals, has_specified_dihedral, dihedral_spec)
            connection = Connection(*connect_vals, has_specified_dihedral=has_spec_dihedral,
                                    specified_dihedral=dihedral_spec)
            self.connections.append(connection)
        # print([x.specified_dihedral for x in self.connections])
        # print([x.has_specified_dihedral for x in self.connections])
        # quit()

    # Read connect file and stores connections and building blocks
    def read_connections(self):
        found_bb = False
        found_connects = False

        with open(Utils.connections_paths, 'r') as f:
            lines = f.readlines()

        # check for keywords:  not used
        # n_bb_kw = lines_have_keyword(lines, "[ BUILDING BLOCKS ]")
        # n_con_kw = lines_have_keyword(lines, "[ CONNECTS ]")
        n_con_w_dihe_kw = lines_have_keyword(lines, "[ ADVANCED DIHEDRAL ]")
        if n_con_w_dihe_kw == 1:
            print("Found Connections with specified dihedral in input.")
            found_connects_w_dih = True
        if n_con_w_dihe_kw > 1:
            print("Keyword [ ADVANCED DIHEDRAL ] appears multiple times in connection file.")
            quit()
        con_lines_w_dihedral = []
        for lidx, line in enumerate(lines):
            if is_comment(line):
                pass
            if "[ BUILDING BLOCKS ]" in line:
                if found_bb:
                    print("Keyword [ BUILDING BLOCKS ] appears multiple times in connection file.")
                    quit()
                found_bb = True
                bb_lines = self.read_building_blocks(lines, lidx)
                self.process_building_blocks(bb_lines)
            if "[ CONNECTS ]" in line:
                if found_connects:
                    print("Keyword [ CONNECTS ] appears multiple times in connection file.")
                    quit()
                found_connects = True
                con_lines = self.read_connects(lines, lidx)
            if "[ ADVANCED DIHEDRAL ]" in line:
                con_lines_w_dihedral = self.read_connects_w_dihedral(lines, lidx)

        self.process_connects(con_lines, con_lines_w_dihedral)

    # read and store the building blocks files paths
    def read_bbs(self):
        for path in Utils.bbs_paths:
            bb_factory = BuildingBlockFactory(path)
            self.bb_factories[bb_factory.resName] = bb_factory

    def connect_bb(self, bb):
        atoms     = copy.copy(bb.get_atom_list())
        bonds     = copy.copy(bb.get_bond_list())
        angles    = copy.copy(bb.get_angle_list())
        dihedrals = copy.copy(bb.get_dihedral_list())

        self.atom_list.extend(atoms)
        self.bond_list.extend(bonds)
        self.angle_list.extend(angles)
        self.dihedral_list.extend(dihedrals)

    def connect_bbs(self):
        for connection in self.connections:
            # print(self.building_blocks)
            firstBB = self.building_blocks[connection.b1]
            secondBB = self.building_blocks[connection.b2]

            atom1 = firstBB.atom_list[connection.atom_b1-1]
            atom2 = secondBB.atom_list[connection.atom_b2-1]
            if connection.has_specified_dihedral:
                specified_dihed = connection.specified_dihedral
                # find partners here for dihedral
                atom0 = firstBB.atom_list[int(specified_dihed["atom_nrs"][0]) - 1]
                atom3 = secondBB.atom_list[int(specified_dihed["atom_nrs"][3]) - 1]
                # pass on dihedral info to all atoms in dihedral
                dihed_atom_nrs = [atom.get_nr() for atom in [atom0, atom1, atom2, atom3]]
                # print(specified_dihed["atom_nrs"])
                # print(connection.atom_b2, atom2.get_nr())
                # print(dihed_atom_nrs)
                for aidx, atom in enumerate([atom0, atom1, atom2, atom3]):
                    atom.set_is_in_specified_dihedral(True)
                    atom.set_specified_dihedral_vals({"idx_in_dihed": aidx, "dihed_atom_nrs": dihed_atom_nrs,
                                                      "dihed_atoms": [atom0, atom1, atom2, atom3],
                                                      "funct": specified_dihed["funct"], "ph0": specified_dihed["ph0"],
                                                      "use_for_lookup": aidx == 1})  # only use atom1 for lookup

                # for atom in [atom0, atom1, atom2, atom3]:
                #     print(atom.get_is_in_specified_dihedral(), atom.get_nr(), atom.get_neighbor_list())

            bond = Bond(atom1, atom2, 2, None)
            self.bond_list.append(bond)


    def __init__(self):
        self.connections  = []
        self.bb_factories = {}
        self.dict_number_resname = {}
        self.atom_list      = []
        self.bond_list      = []
        self.angle_list     = []
        self.dihedral_list  = []
        self.exclusion_list = []
        self.building_blocks = {}
