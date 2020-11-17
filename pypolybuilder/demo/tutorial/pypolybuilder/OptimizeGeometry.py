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

class OptimizeGeometry(object):
    

    def __init__(self, top):
        self.top = top
        self.initForces()
        self.totalForce = 0.0
        self.minNumberOfSteps = 300
        self.numberOfSteps = 600
        self.stepLength = 0.01
        self.maxstepLength = 0.05
        self.minForceDiff = 0.00001
        self.totalBonds = 0.0
        self.totalLJ = 0.0
        self.totalAngles = [0.0,0.0,0.0]
        self.maxForce = 100
        self.currstep = 0
        self.calcPairListParam = 10

    def initForces(self):
        self.forces = []
        vec = [0.0, 0.0, 0.0]
        for atoms in self.top.get_atom_list():
            self.forces.append(vec)

    def getStepToCalcPair(self,multiplier):
        return multiplier * self.calcPairListParam

    def run(self):
        i = 0
        pairMultiplier = 1
        self.initForces()
        self.calcBondedTerms()
        self.calcAngledTerms()
        self.calcLJTerms()
        self.calcTotalForce()
        

        while(i < self.numberOfSteps):
            self.moveAtoms()
            if(self.breakCondition()):
                break
            self.initForces()
            self.calcBondedTerms()
            self.calcAngledTerms()
            self.calcDihedralTerms()
            self.calcLJTerms()
            i = i + 1
            self.currstep = i
            if((i % self.calcPairListParam) == 0):
                self.top.standardPairList()

    #Returns whether the difference of old and new forces absolute values is lower than the minimum allowed
    def breakCondition(self):
        oldTotalForce = self.totalForce
        self.calcTotalForce()

        if (self.currstep > self.minNumberOfSteps):
            if (self.totalForce < oldTotalForce):
                self.stepLength = self.stepLength * 1.2
                if(self.stepLength > (self.maxstepLength)):
                    self.stepLength = self.maxstepLength
            else:
                self.stepLength = self.stepLength * 0.5
        

        totalDiff = oldTotalForce - self.totalForce
        totalDiff = abs(totalDiff)
        if( totalDiff == 0):
            return False
        return (totalDiff < self.minForceDiff)

    #Loop to move all atoms
    def moveAtoms(self):
        atoms = self.top.get_atom_list()
        for atom in atoms:
            self.moveAtom(atom)

    #Change atom coordinates based on the force exerted on it
    def moveAtom(self,atom):
        forceOnAtom = self.forces[int(atom.get_nr())-1]
        if (self.calcVectorModule(forceOnAtom) > self.maxForce):
            forceOnAtom = forceOnAtom * (self.maxForce / self.calcVectorModule(forceOnAtom))
        moveValue = forceOnAtom * self.stepLength
        atom.add_to_x(moveValue[0])
        atom.add_to_y(moveValue[1])
        atom.add_to_z(moveValue[2])

    #Calculate the sum of the forces absolute values
    def calcTotalForce(self):
        self.totalForce = 0.0
        for force in self.forces:
            self.totalForce = self.totalForce + self.calcVectorModule(force)

    def calcBondedForce(self,bond):
        k  = 10000
        b0 = 0.145
        atom1 = bond.get_a_1()
        atom2 = bond.get_a_2()

        if((atom1.get_atomtype() == "H" or atom1.get_atomtype() == "HC") or  (atom2.get_atomtype() == "H" or atom2.get_atomtype() == "HC")):
            k = 15000
            b0 = 0.1

        xdiff = atom1.get_x() - atom2.get_x()
        ydiff = atom1.get_y() - atom2.get_y()
        zdiff = atom1.get_z() - atom2.get_z()
        vec12 = [xdiff,ydiff,zdiff]
        distance = atom1.calcDist(atom2)
        diff = distance - b0
        tempscalar = (-k * diff / distance)
        force = np.multiply(vec12,tempscalar)
        return force

    

    def calcAngledForce(self, angle):
        k = 1000
        cost0 = -0.5  # 120 degree
        atom1 = angle.get_a_1()
        atom2 = angle.get_a_2()
        atom3 = angle.get_a_3()

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

    def calcDihedralForce(self, dihedral):
        if dihedral.get_conformation() is None:
            return

        cosdelta = 1.0
        if dihedral.get_conformation() == 'cis':
            cosdelta = -1.0

        atom1 = dihedral.get_a_1()
        atom2 = dihedral.get_a_2()
        atom3 = dihedral.get_a_3()
        atom4 = dihedral.get_a_4()

        rij = vectorFromAtoms(dihedral.get_a_1(),dihedral.get_a_2())
        rkj = vectorFromAtoms(dihedral.get_a_3(),dihedral.get_a_2())
        rkl = vectorFromAtoms(dihedral.get_a_3(),dihedral.get_a_4())


        k = 33.5

        #rmj = np.cross(rij,rkj)
        #rnk = np.cross(rkj,rkl)

        dkj2 = np.linalg.norm(rkj)*np.linalg.norm(rkj)

        frim = np.dot(rij,rkj) / dkj2
        frln = np.dot(rkl,rkj) / dkj2

        rim = np.array(rij) - frim * np.array(rkj)
        rln = frln * np.array(rkj) - np.array(rkl)

        dim = math.sqrt(np.linalg.norm(rim)*np.linalg.norm(rim))
        dln = math.sqrt(np.linalg.norm(rln)*np.linalg.norm(rln))

        ip = np.dot(rim,rln)


        cosphi = ip / (dim * dln)
        cosphi2 = cosphi*cosphi
        cosphi3 = cosphi2*cosphi
        cosphi4 = cosphi3*cosphi

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



    def calcBondedTerms(self):
        for bond in self.top.get_bond_list():
            atom1 = bond.get_a_1()
            atom2 = bond.get_a_2()
            force = self.calcBondedForce(bond)
            self.forces[int(atom1.get_nr())-1] += force
            self.forces[int(atom2.get_nr())-1] -= force
            self.totalBonds = self.totalBonds + self.calcVectorModule(force)


    def calcAngledTerms(self):
        for angle in self.top.get_angle_list():
            self.calcAngledForce(angle)

    def calcDihedralTerms(self):
        for dihedral in self.top.get_dihedral_list():
            self.calcDihedralForce(dihedral)

    def calcLJForce(self,atom1,atom2):
        c12  = 0.000001421
        c6   = 0.0017
        vecAB = self.vectorFromAtoms(atom1, atom2)
        distance = atom1.calcDist(atom2)
        k = (1/distance) * ((-c6/(distance**6)) + (c12/(distance**12)))
        force = np.multiply(vecAB, k)
        return force

    def calcLJTerms(self):
        for pair in self.top.get_pair_list():
            force = self.calcLJForce(pair[0],pair[1])
            self.totalLJ = self.totalLJ + self.calcVectorModule(force)

            self.forces[int(pair[0].get_nr()) - 1] += force
            self.forces[int(pair[1].get_nr()) - 1] -= force

    def calcVectorModule(self, v):
        return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

    def unit_vector(self,vector):
        return vector / np.linalg.norm(vector)

    def angle_between(self,v1, v2):
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

    def cosine_similarity(self,v1, v2):
        "compute cosine similarity of v1 to v2: (v1 dot v2)/{||v1||*||v2||)"
        sumxx, sumxy, sumyy = 0, 0, 0
        for i in range(len(v1)):
            x = v1[i];
            y = v2[i]
            sumxx += x * x
            sumyy += y * y
            sumxy += x * y
        return sumxy / math.sqrt(sumxx * sumyy)

    def cosine_2unit_v(self, v1,v2):
        scalar = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
        return scalar


    def vectorFromAtoms(self,atom1,atom2):
        vec12 = [atom1.get_x() - atom2.get_x(), atom1.get_y() - atom2.get_y(),
                 atom1.get_z() - atom2.get_z()]
        return vec12