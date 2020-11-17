#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 12:13:44 2018

@author: vitor
"""

import copy
from random import uniform
# Doug
import numpy as np
from random import uniform, random, randint


class VBGA_U:
    def __init__(self, individualClass, fitnessMethod, selectionMethod, crossoverMethod, crossoverRate, mutationMethod,
                 mutationRate, populationSize):
        self.populationSize = populationSize
        self.individualClass = individualClass
        self.fitnessMethod = fitnessMethod
        self.selectionMethod = selectionMethod
        self.crossoverMethod = crossoverMethod
        self.crossoverRate = crossoverRate
        self.mutationRate = mutationRate
        self.mutationMethod = mutationMethod
        self.population = []

    def createRandomIndividual(self):
        ind = self.individualClass()
        ind.randomize()
        return ind

    def createInitialPopulation(self):
        self.population = []
        for i in range(0, self.populationSize):
            ind = self.createRandomIndividual()
            self.population.append(ind)

    def evaluateIndividual(self, individual):
        individual.fitValue = self.fitnessMethod(individual)

    def evaluatePopulation(self):
        for individual in self.population:
            self.evaluateIndividual(individual)

    def averageFitness(self):
        max = 0
        for individual in self.population:
            # print(individual.fitValue)
            if individual.fitValue > max:
                max = individual.fitValue
        return max

    # Luigi
    def normalizationFitness(self):
        sumFitness = 0
        for individual in self.population:
                sumFitness += individual.fitValue
        for individual in self.population:
            individual.fitValue = individual.fitValue / sumFitness

    def applySelectionMethod(self):
        selectedPopulation = self.selectionMethod(self.population)
        return selectedPopulation

    # Luigi
    def elitism(self, selectedPopulation):
        nextPopulation = copy.deepcopy(selectedPopulation)
        self.population.sort(key=lambda x: x.fitValue, reverse=True)
        for i in range(0, 4):
            ind = self.population[i]
            nextPopulation.append(ind)

        self.population = nextPopulation

    def fillPopulation(self, selectedPopulation):
        nextPopulation = copy.deepcopy(selectedPopulation)
        for i in range(len(selectedPopulation), self.populationSize):
            ind = self.createRandomIndividual()
            nextPopulation.append(ind)

        self.population = nextPopulation

    # Luigi
    def applyCrossoverMethod(self, selected):
        # populationCross = copy.deepcopy(selected)
        populationCross = []
        for i in range(0, len(selected)):
            for j in range(i + 1, len(selected)):
                dice = uniform(0, 100)
                if dice < self.crossoverRate:
                    father = selected[i]
                    mother = selected[j]
                    childOne = self.individualClass()
                    childOne.randomize()
                    childTwo = self.individualClass()
                    childTwo.randomize()
                    children = self.crossoverMethod(father, mother, childOne, childTwo)
                    # children[0] = self.applyMutationMethod(children[0])
                    # children[1] = self.applyMutationMethod(children[1])
                    # populationCross = populationCross + children  # Elitism - fathers + children
                    populationCross.append(children[0])
                    populationCross.append(children[1])
        # print("DEBUG -- population size: ", len(populationCross))
        return populationCross

    # Doug
    def applyCrossoverMethodTour(self, selected):
        populationCross = copy.deepcopy(selected)
        dice = uniform(0, 100)
        if dice < self.crossoverRate:
            father = selected[0]
            mother = selected[1]
            childOne = self.individualClass()
            childTwo = self.individualClass()
            children = self.crossoverMethod(father, mother, childOne, childTwo)
            children[0] = self.applyMutationMethod(children[0])
            children[1] = self.applyMutationMethod(children[1])
            populationCross = [childOne, childTwo]
        return populationCross

    def applyMutationMethod(self, individual):
        dice = uniform(0, 100)
        if dice < self.mutationRate:
            individual = self.mutationMethod(individual)
        return individual

    # Luigi
    def mutationPopulation(self, selectedPopulation):
        for individual in selectedPopulation:
            self.applyMutationMethod(individual)

    def getBest(self):
        max = 0
        best = None
        for individual in self.population:
            # print(individual.fitValue)
            if individual.fitValue > max:
                best = individual
                max = individual.fitValue
        # Doug
        print("Best fitness: ", max)
        return best

    def run(self, maxIterations=100):
        self.createInitialPopulation()
        self.evaluatePopulation()  # Update fitness
        for i in range(0, maxIterations):
            self.normalizationFitness()
            selected = self.applySelectionMethod()  # Selection of the fathers
            populationCross = self.applyCrossoverMethod(selected)  # Create childen - populationCross is formed for fathers + chlidren
            self.mutationPopulation(populationCross)
            self.elitism(populationCross)  # Best 4 individual in next populition
            self.fillPopulation(self.population)  # Complite populition with randomzide individuals
            # self.fillPopulation(selected)
            self.evaluatePopulation()  # Update fitness
            avFitness = self.averageFitness()  # Best fitness
            print(avFitness)

    # Doug
    def runTour(self, maxIterations=100):
        self.createInitialPopulation()
        for i in range(0, maxIterations):
            a = np.random.randint(self.populationSize)
            b = np.random.randint(self.populationSize)
            c = np.random.randint(self.populationSize)
            d = np.random.randint(self.populationSize)
            individual_a = self.population[a]
            individual_b = self.population[b]
            individual_c = self.population[c]
            individual_d = self.population[d]
            self.evaluateIndividual(individual_a)
            self.evaluateIndividual(individual_b)
            self.evaluateIndividual(individual_c)
            self.evaluateIndividual(individual_d)
            winner_ab = a
            looser_ab = b
            if individual_a.fitValue < individual_b.fitValue:
                winner_ab = b
                looser_ab = a
            winner_cd = c
            looser_cd = d
            if individual_c.fitValue < individual_d.fitValue:
                winner_cd = d
                looser_cd = c
            selected = [self.population[winner_ab], self.population[winner_cd]]
            populationCross = self.applyCrossoverMethodTour(selected)
            self.population[looser_ab] = copy.deepcopy(populationCross[0])
            self.population[looser_cd] = copy.deepcopy(populationCross[1])


def DefaultSelectionMethod(population, selectionSize=10):
    population.sort(key=lambda x: x.fitValue, reverse=True)
    return population[:selectionSize]  # return 0,1,2....selectionSize


# Luigi
def RoulettSelectionMethod(population, selectionSize=10):  # FitValue must be normalized
    select = []
    while len(select) < selectionSize:
        roulettPoint = random()  # (0.0, 1.0)
        fitPoint = 0
        for individual in population:
            fitPoint += individual.fitValue  # FitPoint should be incremented until it is larger than roulettPoint
            if roulettPoint < fitPoint:
                select.append(individual)
                break
    return select  # return 0,1,2....selectionSize


def DefaultCrossoverMethod(father, mother, childOne, childTwo):
    child.x = (father.x + mother.x) / 2
    child.y = (father.y + mother.y) / 2

    return child


def UniformCrossoverMethod(father, mother, childOne, childTwo):
    for i in range(0, len(father.matrix.dihvalue)):
        childOne.matrix.dihvalue[i] = (father.matrix.dihvalue[i] + mother.matrix.dihvalue[
            i]) / 2  # I changed the 2 for 2.5 - Luigi
        childTwo.matrix.dihvalue[i] = (father.matrix.dihvalue[i] + mother.matrix.dihvalue[
            i]) / 3  # I changed the 3 for 1.5 - Luigi
    return [childOne, childTwo]


# Jorge
def HeuristicCrossoverMethod(father, mother, childOne, childTwo):
    r = random()
    if father.fitValue >= mother.fitValue:
        for i in range(0, len(father.matrix.dihvalue)):
            childOne.matrix.dihvalue[i] = father.matrix.dihvalue[i] + r * (
                        father.matrix.dihvalue[i] - mother.matrix.dihvalue[i])
            childTwo.matrix.dihvalue[i] = father.matrix.dihvalue[i]
    else:
        for i in range(0, len(father.matrix.dihvalue)):
            childOne.matrix.dihvalue[i] = r * (mother.matrix.dihvalue[i] - father.matrix.dihvalue[i]) + \
                                          mother.matrix.dihvalue[i]
            childTwo.matrix.dihvalue[i] = mother.matrix.dihvalue[i]

    return [childOne, childTwo]


# Jorge
def AritmeticCrossoverMethod(father, mother, childOne, childTwo):
    r = random()
    for i in range(0, len(father.matrix.dihvalue)):
        childOne.matrix.dihvalue[i] = r * father.matrix.dihvalue[i] + (1 - r) * mother.matrix.dihvalue[i]
        childTwo.matrix.dihvalue[i] = (1 - r) * father.matrix.dihvalue[i] + r * mother.matrix.dihvalue[i]
    return [childOne, childTwo]


# Luigi
def MatingCrossoverMethod(father, mother, childOne, childTwo, method='Two Points'):
    if method == 'Single Point':
        pivot_point = randint(1, len(father.matrix.dihvalue))
        childOne.matrix.dihvalue = father.matrix.dihvalue[0:pivot_point] + mother.matrix.dihvalue[pivot_point:]
        childTwo.matrix.dihvalue = mother.matrix.dihvalue[0:pivot_point] + father.matrix.dihvalue[pivot_point:]
    if method == 'Two Points':
        pivot_point_1 = randint(1, len(father.matrix.dihvalue) - 1)
        pivot_point_2 = randint(pivot_point_1 + 1, len(father.matrix.dihvalue))
        childOne.matrix.dihvalue = father.matrix.dihvalue[0:pivot_point_1] + mother.matrix.dihvalue[pivot_point_1:pivot_point_2] + father.matrix.dihvalue[pivot_point_2:]
        childTwo.matrix.dihvalue = mother.matrix.dihvalue[0:pivot_point_1] + father.matrix.dihvalue[pivot_point_1:pivot_point_2] + mother.matrix.dihvalue[pivot_point_2:]
    return [childOne, childTwo]


def DefaultMutation(individual):
    probability = 0.1
    for i in range(0, len(individual.matrix.dihvalue)):
        dice = uniform(0, 1)
        # print("DEBUG -- prob and dice ", probability, " ", dice)
        if dice < probability:
            a = np.random.randint(12)
            newangle = float(a * 30)
            individual.matrix.dihvalue[i] = newangle
    return individual


def calcRgy(individual):
    natoms = len(individual.matrix.atoms)
    xyzarr = individual.matrix.converToXYZ()
    x, y, z = 0.0, 0.0, 0.0
    for i in range(natoms):
        x += xyzarr[i][0]
        y += xyzarr[i][1]
        z += xyzarr[i][2]
        # print ('{:<4s}\t{:>11.5f}\t{:>11.5f}\t{:>11.5f}'.format(self.atoms[i].get_atomtype()[0], xyzarr[i][0], xyzarr[i][1], xyzarr[i][2]))
    centerOfGeom = np.array([x, y, z]) / float(natoms)
    rcum = 0.0
    for i in range(natoms):
        rcum += np.linalg.norm(xyzarr[i] - centerOfGeom)
    return (rcum / natoms)


def calcSphericity(individual):
    natoms = len(individual.matrix.atoms)
    xyzarr = individual.matrix.converToXYZ()
    w = [1/3, 1/3, 1/3]  # Weights, define final shape
    x, y, z = 0.0, 0.0, 0.0
    for i in range(natoms):
        x += xyzarr[i][0]
        y += xyzarr[i][1]
        z += xyzarr[i][2]
        # print ('{:<4s}\t{:>11.5f}\t{:>11.5f}\t{:>11.5f}'.format(self.atoms[i].get_atomtype()[0], xyzarr[i][0], xyzarr[i][1], xyzarr[i][2]))
    centerOfGeom = np.array([x, y, z]) / float(natoms)  # Centroid
    rcum = 0.0
    rVec = [0.0, 0.0, 0.0]
    for i in range(natoms):
        rcum += np.linalg.norm(xyzarr[i] - centerOfGeom)  # Euclidean norm
        rVec += (xyzarr[i][0] - centerOfGeom)
    sphericity = w[0] * rVec[0] / rcum + w[1] * rVec[1] / rcum + w[2] * rVec[2] / rcum
    return sphericity


def calcDistSum(individual):
    natoms = len(individual.matrix.atoms)
    xyzarr = individual.matrix.converToXYZ()
    x, y, z = 0.0, 0.0, 0.0
    distSum = 0.0
    for i in range(natoms):
        x_i = xyzarr[i][0]
        y_i = xyzarr[i][1]
        z_i = xyzarr[i][2]
        for j in range(i, natoms):
            x_j = xyzarr[j][0]
            y_j = xyzarr[j][1]
            z_j = xyzarr[j][2]
            dist = (x_i - x_j) + (y_i - y_j) + (z_i - z_j)
            dist = dist * dist
            distSum += dist
    return distSum