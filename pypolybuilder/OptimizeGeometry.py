'''
Created on Apr 20, 2015

@author: vitor
'''

from __future__ import division
import sys
import time
import numpy as np
import math
from Polymer import Polymer
from Dendrimer import Dendrimer
from Utils import unit_vector
from Utils import vectorFromAtoms
from Utils import nsteps, nskipLJ, stepLength
from Utils import get_ff_params, forcefield_path
import copy
import random

sys.setrecursionlimit(10000)


class OptimizeGeometry(object):

    def __init__(self, top):
        self.top = top
        self.bestTop = top
        self.initForces()
        self.initOldMoveValue()
        self.initStepLengthVec()
        self.oldForces = []
        self.initOldPosition()
        self.totalForce = 0.0
        self.bestForce = 0
        self.minNumberOfSteps = 20
        self.maxNumberOfSteps = 600
        self.numberOfSteps = nsteps
        self.nskipLJ = nskipLJ
        self.stepLength = stepLength
        self.maxstepLength = 0.05
        self.minForceDiff = 0.00001
        self.totalBonds = 0.0
        self.totalLJ = 0.0
        self.totalAngles = [0.0, 0.0, 0.0]
        self.maxForce = 100
        self.currstep = 0
        self.calcPairListParam = 10
        self.referenceAngle = 0
        self.technical = 0

    def initOldPosition(self):
        self.oldPosition = []
        vec = [0.0, 0.0, 0.0]
        for atoms in self.top.get_atom_list():
            self.oldPosition.append(vec)

    def initForces(self):
        self.forces = []
        vec = [0.0, 0.0, 0.0]
        for atoms in self.top.get_atom_list():
            self.forces.append(vec)

    def initOldMoveValue(self):
        self.oldMoveValue = []
        vec = [0.0, 0.0, 0.0]
        for atoms in self.top.get_atom_list():
            self.oldMoveValue.append(vec)

    def initStepLengthVec(self):
        self.stepLengthVec = []
        vec = [0.00001, 0.00001, 0.00001]
        for atoms in self.top.get_atom_list():
            self.stepLengthVec.append(vec)

    def getStepToCalcPair(self, multiplier):
        return multiplier * self.calcPairListParam

    def printAbsForces(self, i):
        print(
            f'{i:<10}{self.totalBonds:<25}{self.totalAngles[0]:<25}{self.totalAngles[1]:<25}{self.totalAngles[2]:<25}'
            f'{self.totalLJ:<25}{self.totalForce:<25}')
        self.totalBonds = 0.0
        self.totalLJ = 0.0
        self.totalAngles = [0.0, 0.0, 0.0]

    def run(self):

        ff_bonds, ff_angles, ff_dihedrals, ff_LJs = get_ff_params(forcefield_path)

        print(
            f'{"Step":<10}{"totalBonds":<25}{"totalAngles[0]":<25}{"totalAngles[1]":<25}{"totalAngles[2]":<25}'
            f'{"totalLJ":<25}{"totalForce":<25}')
        i = 0
        global totalF
        global totalE
        totalF = []
        totalE = []
        pairMultiplier = 1
        self.initForces()
        self.initOldPosition()
        self.calcBondedTerms(ff_bonds)
        self.calcAngledTerms(ff_angles)
        if self.nskipLJ == 0:
            self.calcLJTerms(ff_LJs)
        self.calcTotalForce()
        totalEnergy = self.calcTotalEnery(ff_bonds, ff_angles, ff_dihedrals, ff_LJs)
        totalF.append(self.totalForce)
        totalE.append(totalEnergy)
        self.referenceAngle = self.totalForce / self.numberOfSteps
        self.bestForce = self.totalForce

        while i < self.numberOfSteps:
            self.moveAtoms()
            if self.breakCondition(ff_bonds, ff_angles, ff_dihedrals, ff_LJs):
                break
            self.initForces()
            self.calcBondedTerms(ff_bonds)
            self.calcAngledTerms(ff_angles)
            self.calcDihedralTerms(ff_dihedrals)
            if i > self.nskipLJ:
                self.calcLJTerms(ff_LJs)
            if i%(0.1*self.numberOfSteps) == 0:
                self.printAbsForces(i)
            i = i + 1
            self.currstep = i
            if i % self.calcPairListParam == 0:
                self.top.standardPairList()

        self.top = copy.deepcopy(self.bestTop)
        print('\nBest totalForce: {}'.format(self.bestForce))

        # totalf_df = pd.DataFrame(totalF)
        # totalf_df.to_csv('totalForce.csv')

        # totalE_df = pd.DataFrame(totalE)
        # totalE_df.to_csv('totalEnergy.csv')

    # Returns whether the difference of old and new forces absolute values is lower than the minimum allowed
    # @property
    def breakCondition(self, ff_bonds, ff_angles, ff_dihedrals, ff_LJs):
        oldTotalForce = self.totalForce
        self.oldForces = copy.deepcopy(self.forces)
        self.calcTotalForce()
        totalEnergy = self.calcTotalEnery(ff_bonds, ff_angles, ff_dihedrals, ff_LJs)

        newAngle = oldTotalForce - self.totalForce

        totalF.append(self.totalForce)
        totalE.append(totalEnergy)

        if self.totalForce < self.bestForce:
            self.bestForce = self.totalForce
            self.bestTop = copy.deepcopy(self.top)

        if self.referenceAngle > newAngle > 0:  # Counts the number of worsening descents
            self.technical += 1

        if self.technical >= 15:
            if self.totalForce < oldTotalForce:
                self.stepLength = self.stepLength * 1.4
                if self.stepLength > self.maxstepLength:
                    self.stepLength = self.maxstepLength
            else:
                self.stepLength = self.stepLength * 0.5

        totalDiff = oldTotalForce - self.totalForce
        totalDiff = abs(totalDiff)
        if totalDiff == 0:
            return False
        return totalDiff < self.minForceDiff

    # Loop to move all atoms
    def moveAtoms(self):
        atoms = self.top.get_atom_list()
        for atom in atoms:
            self.moveAtom(atom)

    # Change atom coordinates based on the force exerted on it
    def moveAtom(self, atom):
        forceOnAtom = self.forces[int(atom.get_nr()) - 1]
        self.oldPosition = [atom._Atom__x, atom._Atom__y, atom._Atom__z]
        if (self.calcVectorModule(forceOnAtom) > self.maxForce):
            forceOnAtom = forceOnAtom * (self.maxForce / self.calcVectorModule(forceOnAtom))
        moveValue = forceOnAtom * self.stepLength
        atom.add_to_x(moveValue[0])
        atom.add_to_y(moveValue[1])
        atom.add_to_z(moveValue[2])

    def moveAtomsReturn(self):
        atoms = self.top.get_atom_list()
        for atom in atoms:
            self.moveAtom(atom)

    # Change atom coordinates based on the force exerted on it
    def moveAtomReturn(self, atom):
        atom.set_x(self.oldPosition[0])
        atom.set_y(self.oldPosition[1])
        atom.set_z(self.oldPosition[2])

    # Calculate the sum of the forces absolute values
    def calcTotalForce(self):
        self.totalForce = 0.0
        for force in self.forces:
            self.totalForce = self.totalForce + self.calcVectorModule(force)

    def calcBondedForce(self, bond, ff_bonds):
        k = 10000
        b0 = 0.145
        
        atom1 = bond.get_a_1()
        atom2 = bond.get_a_2()

        if ((atom1.get_atomtype() == "H" or atom1.get_atomtype() == "HC") or (
                atom2.get_atomtype() == "H" or atom2.get_atomtype() == "HC")):
            k = 15000
            b0 = 0.1

        # Overwrite the hardcoded params
        if ff_bonds is not None:
            b0 = ff_bonds[bond.get_param()][0]
            k  = ff_bonds[bond.get_param()][1]

        xdiff = atom1.get_x() - atom2.get_x()
        ydiff = atom1.get_y() - atom2.get_y()
        zdiff = atom1.get_z() - atom2.get_z()
        vec12 = [xdiff, ydiff, zdiff]
        distance = atom1.calcDist(atom2)
        diff = distance - b0
        tempscalar = (-k * diff / distance)
        force = np.multiply(vec12, tempscalar)
        return force

    def calcAngledForce(self, angle, ff_angles):
        k = 1000
        cost0 = -0.5  # 120 degree

        atom1 = angle.get_a_1()
        atom2 = angle.get_a_2()
        atom3 = angle.get_a_3()

        # Overwrite the hardcoded params
        if ff_angles is not None:
            t0    = ff_angles[angle.get_param()][0]
            cost0 = np.cos(t0*(np.pi/180))
            # k     = ff_angles[angle.get_param()][1]

        vecAB = self.vectorFromAtoms(atom1, atom2)
        vecAB = unit_vector(vecAB)
        vecCB = vectorFromAtoms(atom3, atom2)
        vecCB = unit_vector(vecCB)
        costeta = self.cosine_2unit_v(vecAB, vecCB)
        df = -k * (costeta - cost0)
        forceA = df * (vecCB - vecAB * costeta)
        forceC = df * (vecAB - vecCB * costeta)
        forceB = -1 * forceA - forceC

        self.forces[int(atom1.get_nr()) - 1] += forceA
        self.forces[int(atom2.get_nr()) - 1] += forceB
        self.forces[int(atom3.get_nr()) - 1] += forceC
        self.totalAngles[0] += self.calcVectorModule(forceA)
        self.totalAngles[1] += self.calcVectorModule(forceB)
        self.totalAngles[2] += self.calcVectorModule(forceC)

    def calcDihedralForce(self, dihedral, ff_dihedrals):
        if dihedral.get_conformation() is None:
            return

        cosdelta = 1.0
        if dihedral.get_conformation() == 'cis':
            cosdelta = -1.0

        atom1 = dihedral.get_a_1()
        atom2 = dihedral.get_a_2()
        atom3 = dihedral.get_a_3()
        atom4 = dihedral.get_a_4()

        # Overwrite the hardcoded params
        if ff_dihedrals is not None:
            delta    = ff_dihedrals[dihedral.get_param()][0]
            cosdelta = np.cos(delta*(np.pi/180))

        rij = vectorFromAtoms(dihedral.get_a_1(), dihedral.get_a_2())
        rkj = vectorFromAtoms(dihedral.get_a_3(), dihedral.get_a_2())
        rkl = vectorFromAtoms(dihedral.get_a_3(), dihedral.get_a_4())

        k = 33.5
        # k = 100000

        # rmj = np.cross(rij,rkj)
        # rnk = np.cross(rkj,rkl)

        dkj2 = np.linalg.norm(rkj) * np.linalg.norm(rkj)

        frim = np.dot(rij, rkj) / dkj2
        frln = np.dot(rkl, rkj) / dkj2

        rim = np.array(rij) - frim * np.array(rkj)
        rln = frln * np.array(rkj) - np.array(rkl)

        dim = math.sqrt(np.linalg.norm(rim) * np.linalg.norm(rim))
        dln = math.sqrt(np.linalg.norm(rln) * np.linalg.norm(rln))

        ip = np.dot(rim, rln)

        cosphi = ip / (dim * dln)
        cosphi2 = cosphi * cosphi
        cosphi3 = cosphi2 * cosphi
        cosphi4 = cosphi3 * cosphi

        cosmphi = cosphi
        dcosmphi = 1

        ki = -k * cosdelta * dcosmphi / dim
        kl = -k * cosdelta * dcosmphi / dln
        kj1 = frim - 1.0
        kj2 = frln

        fi = ki * (np.array(rln) / dln - np.array(rim) / dim * cosphi)
        fl = kl * (np.array(rim) / dim - np.array(rln) / dln * cosphi)
        fj = kj1 * np.array(fi) - kj2 * np.array(fl)
        fk = -1.0 * (np.array(fi) + np.array(fl) + np.array(fj))

        self.forces[int(atom1.get_nr()) - 1] += fi
        self.forces[int(atom2.get_nr()) - 1] += fj
        self.forces[int(atom3.get_nr()) - 1] += fk
        self.forces[int(atom4.get_nr()) - 1] += fl

    def calcBondedTerms(self, ff_bonds):
        for bond in self.top.get_bond_list():
            atom1 = bond.get_a_1()
            atom2 = bond.get_a_2()
            force = self.calcBondedForce(bond, ff_bonds)
            self.forces[int(atom1.get_nr()) - 1] += force
            self.forces[int(atom2.get_nr()) - 1] -= force
            self.totalBonds = self.totalBonds + self.calcVectorModule(force)

    def calcAngledTerms(self, ff_angles):
        for angle in self.top.get_angle_list():
            self.calcAngledForce(angle, ff_angles)

    def calcDihedralTerms(self, ff_dihedrals):
        for dihedral in self.top.get_dihedral_list():
            self.calcDihedralForce(dihedral, ff_dihedrals)

    def calcLJForce(self, pair, ff_LJs):
        c12 = 0.000001421
        c6 = 0.0017

        # Overwrite the hardcoded params
        if ff_LJs is not None:
            atom1_name = pair[0].get_atomtype()
            atom2_name = pair[1].get_atomtype()
            atoms=f'{atom1_name} {atom2_name}'
            if atoms in ff_LJs:
                c6  = ff_LJs[atoms][0]
                c12 = ff_LJs[atoms][1]

        vecAB = self.vectorFromAtoms(pair[0], pair[1])
        distance = pair[0].calcDist(pair[1])
        # k = (1/distance) * ((-c6/(distance**6)) + (c12/(distance**12)))
        k = (12 * c12 / (distance ** 11))
        force = np.multiply(vecAB, k)
        return force

    def calcLJTerms(self, ff_LJs):
        for pair in self.top.get_pair_list():
            # print str(pair[0].get_nr()) + "   " + str(pair[1].get_nr())
            force = self.calcLJForce(pair, ff_LJs)
            self.totalLJ = self.totalLJ + self.calcVectorModule(force)

            self.forces[int(pair[0].get_nr()) - 1] += force
            self.forces[int(pair[1].get_nr()) - 1] -= force

    def calcVectorModule(self, v):
        return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

    def unit_vector(self, vector):
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    def calcMediaAngle(self):
        media = 0.0
        for angle in self.top.get_angle_list():
            atom1 = angle.get_a_1()
            atom2 = angle.get_a_2()
            atom3 = angle.get_a_3()
            vecAB = vectorFromAtoms(atom2, atom1)
            vecAB = unit_vector(vecAB)
            vecBC = vectorFromAtoms(atom2, atom3)
            vecBC = unit_vector(vecBC)
            teta = self.cosine_similarity(vecAB, vecBC)

            media += teta

    def calcMediaTeta(self):
        media = 0.0
        for angle in self.top.get_angle_list():
            atom1 = angle.get_a_1()
            atom2 = angle.get_a_2()
            atom3 = angle.get_a_3()
            vecAB = vectorFromAtoms(atom2, atom1)
            vecAB = unit_vector(vecAB)
            vecBC = vectorFromAtoms(atom2, atom3)
            vecBC = unit_vector(vecBC)
            teta = self.angle_between(vecAB, vecBC)

            media += teta

    def cosine_similarity(self, v1, v2):
        "compute cosine similarity of v1 to v2: (v1 dot v2)/{||v1||*||v2||)"
        sumxx, sumxy, sumyy = 0, 0, 0
        for i in range(len(v1)):
            x = v1[i];
            y = v2[i]
            sumxx += x * x
            sumyy += y * y
            sumxy += x * y
        return sumxy / math.sqrt(sumxx * sumyy)

    def cosine_2unit_v(self, v1, v2):
        scalar = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
        return scalar

    def vectorFromAtoms(self, atom1, atom2):
        vec12 = [atom1.get_x() - atom2.get_x(), atom1.get_y() - atom2.get_y(),
                 atom1.get_z() - atom2.get_z()]
        return vec12

    def calcBondedEnergy(self, bond, ff_bonds):
        k = 10000
        b0 = 0.145
        atom1 = bond.get_a_1()
        atom2 = bond.get_a_2()

        if ((atom1.get_atomtype() == "H" or atom1.get_atomtype() == "HC") or (
                atom2.get_atomtype() == "H" or atom2.get_atomtype() == "HC")):
            k = 15000
            b0 = 0.1

        #Overwrite the hardcoded params
        if ff_bonds is not None:
            b0 = ff_bonds[bond.get_param()][0]
            k  = ff_bonds[bond.get_param()][1]

        b = atom1.calcDist(atom2)
        energy = 0.5 * k * (b - b0) ** 2
        return energy

    def calcAngleEnergy(self, angle, ff_angles):
        k = 1000
        cost0 = -0.5  # 120 degree
        atom1 = angle.get_a_1()
        atom2 = angle.get_a_2()
        atom3 = angle.get_a_3()

        #Overwrite the hardcoded params
        if ff_angles is not None:
            t0    = ff_angles[angle.get_param()][0]
            cost0 = np.cos(t0*(np.pi/180))
            k     = ff_angles[angle.get_param()][1]

        vecAB = self.vectorFromAtoms(atom1, atom2)
        vecAB = unit_vector(vecAB)
        vecCB = vectorFromAtoms(atom3, atom2)
        vecCB = unit_vector(vecCB)
        costeta = self.cosine_2unit_v(vecAB, vecCB)
        energy = 0.5 * k * (costeta - cost0) ** 2
        return energy

    def calcDihedralEnergy(self, dihedral, ff_dihedrals):
        if dihedral.get_conformation() is None:
            return 0

        cosdelta = 1.0
        if dihedral.get_conformation() == 'cis':
            cosdelta = -1.0

        rij = vectorFromAtoms(dihedral.get_a_1(), dihedral.get_a_2())
        rkj = vectorFromAtoms(dihedral.get_a_3(), dihedral.get_a_2())
        rkl = vectorFromAtoms(dihedral.get_a_3(), dihedral.get_a_4())

        # Overwrite the hardcoded params
        if ff_dihedrals is not None:
            delta    = ff_dihedrals[dihedral.get_param()][0]
            cosdelta = np.cos(delta*(np.pi/180))

        k = 33.5

        dkj2 = np.linalg.norm(rkj) * np.linalg.norm(rkj)

        frim = np.dot(rij, rkj) / dkj2
        frln = np.dot(rkl, rkj) / dkj2

        rim = np.array(rij) - frim * np.array(rkj)
        rln = frln * np.array(rkj) - np.array(rkl)

        dim = math.sqrt(np.linalg.norm(rim) * np.linalg.norm(rim))
        dln = math.sqrt(np.linalg.norm(rln) * np.linalg.norm(rln))

        ip = np.dot(rim, rln)

        cosphi = ip / (dim * dln)
        # This is a simplified formula considering only multiplicity 2

        cosmphi = cosphi

        energy = k * (1 + cosdelta * cosmphi)

        return energy

    def calcLJEnergy(self, pair, ff_LJs):
        c12 = 0.000001421
        c6 = 0.0017

        # Overwrite the hardcoded params
        if ff_LJs is not None:
            atom1_name = pair[0].get_atomtype()
            atom2_name = pair[1].get_atomtype()
            atoms=f'{atom1_name} {atom2_name}'
            if atoms in ff_LJs:
                c6  = ff_LJs[atoms][0]
                c12 = ff_LJs[atoms][1]

        distance = pair[0].calcDist(pair[1])
        # k = (1/distance) * ((-c6/(distance**6)) + (c12/(distance**12)))
        energy = c12 / (distance ** 12)
        return energy

    def calcTotalEnery(self, ff_bonds, ff_angles, ff_dihedrals, ff_LJs):
        bondEnergy = 0.0
        angleEnergy = 0.0
        dihedralEnergy = 0.0
        LJEnergy = 0.0

        for bond in self.top.get_bond_list():
            bondEnergy += self.calcBondedEnergy(bond, ff_bonds)

        for angle in self.top.get_angle_list():
            angleEnergy += self.calcAngleEnergy(angle, ff_angles)

        for dihedral in self.top.get_dihedral_list():
            dihedralEnergy += self.calcDihedralEnergy(dihedral, ff_dihedrals)

        for pair in self.top.get_pair_list():
            LJEnergy += self.calcLJEnergy(pair, ff_LJs)

        total = bondEnergy + angleEnergy + dihedralEnergy + LJEnergy

        return total
