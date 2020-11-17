#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 12:13:44 2018

@author: vitor
"""

import copy
from random import uniform
import numpy as np
import pandas as pd
from random import uniform, random, randint
import matplotlib.pyplot as plt


class VBGA:
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
        self.fitValueMatrix = []

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
        self.population.sort(key=lambda x: x.fitValue, reverse=True)
        best = self.population[0].fitValue
        worse = self.population[-1].fitValue
        av = self.averageFitness()
        self.fitValueMatrix.append([best, av, worse])

    def averageFitness(self):
        s = 0
        for individual in self.population:
            s += individual.fitValue
        return s / len(self.population)

    def bestFitness(self):
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

    def applyCrossoverMethod(self, selected):
        populationCross = []
        for i in range(0, len(selected)):
            dice = uniform(0, 100)
            if dice < self.crossoverRate:
                j = randint(0, len(selected) - 1)
                father = selected[i]
                mother = selected[j]
                childOne = self.individualClass()
                childOne.randomize()
                childTwo = self.individualClass()
                childTwo.randomize()
                children = self.crossoverMethod(father, mother, childOne, childTwo)
                populationCross.append(children[0])
                populationCross.append(children[1])
        return populationCross

    def applyMutationMethod(self, individual):
        dice = uniform(0, 100)
        if dice < self.mutationRate:
            individual = self.mutationMethod(individual)
        return individual

    def mutationPopulation(self, selectedPopulation):
        for individual in selectedPopulation:
            self.applyMutationMethod(individual)

    def getBest(self):
        max = 0
        best = None
        for individual in self.population:
            if individual.fitValue > max:
                best = individual
                max = individual.fitValue
        print("Best fitness: ", max)
        return best

    def run(self, maxIterations=100):
        print(f'{"Generation":<15}{"Best Fitness":<25}')
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
            bestFitness = self.bestFitness()  # Best fitness
            print(f'{i:<15}{bestFitness:<25}')
        '''dfFitness = pd.DataFrame({'Bests': np.array(self.fitValueMatrix)[:, 0],
                                  'Averages': np.array(self.fitValueMatrix)[:, 1],
                                  'Worses': np.array(self.fitValueMatrix)[:, 2]})
        dfFitness.plot.line()
        plt.show()
        plt.savefig('agLines.png')'''


def DefaultSelectionMethod(population, selectionSize=10):
    population.sort(key=lambda x: x.fitValue, reverse=True)
    return population[:selectionSize]  # return 0,1,2....selectionSize


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


def UniformCrossoverMethod(father, mother, childOne, childTwo):
    for i in range(0, len(father.matrix.dihvalue)):
        childOne.matrix.dihvalue[i] = (father.matrix.dihvalue[i] + mother.matrix.dihvalue[
            i]) / 2
        childTwo.matrix.dihvalue[i] = (father.matrix.dihvalue[i] + mother.matrix.dihvalue[
            i]) / 3
    return [childOne, childTwo]


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


def AritmeticCrossoverMethod(father, mother, childOne, childTwo):
    r = random()
    for i in range(0, len(father.matrix.dihvalue)):
        if mother.imp[i]:  # If mother.imp[i] is 0, it's improper
            childOne.matrix.dihvalue[i] = r * father.matrix.dihvalue[i] + (1 - r) * mother.matrix.dihvalue[i]
            childTwo.matrix.dihvalue[i] = (1 - r) * father.matrix.dihvalue[i] + r * mother.matrix.dihvalue[i]
    return [childOne, childTwo]


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
        if dice < probability:
            a = np.random.randint(0, 360)
            newangle = float(a)
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
    centerOfGeom = np.array([x, y, z]) / float(natoms)
    rcum = 0.0
    for i in range(natoms):
        rcum += np.linalg.norm(xyzarr[i] - centerOfGeom)
    return (rcum / natoms)


def calcSphericity(individual):
    natoms = len(individual.matrix.atoms)
    xyzarr = individual.matrix.converToXYZ()
    x, y, z = 0.0, 0.0, 0.0
    for i in range(natoms):
        x += xyzarr[i][0]
        y += xyzarr[i][1]
        z += xyzarr[i][2]
    centerOfGeom = np.array([x, y, z]) / float(natoms)

    df = pd.DataFrame({'x': xyzarr[:, 0], 'y': xyzarr[:, 1], 'z': xyzarr[:, 2]})
    rX = 2 * df['x'].std()
    rY = 2 * df['y'].std()
    rZ = 2 * df['z'].std()

    w = [1, 1, 1]

    sphericity = ((rX * rY * rZ) ** (1 / 3)) / max([rX / w[0], rY / w[1], rZ / w[2]]) * 10
    rcum = 0.0
    for i in range(natoms):
        rcum += np.linalg.norm(xyzarr[i] - centerOfGeom)
    return (rcum / natoms) * sphericity


def calcDistSum(individual):
    natoms = len(individual.matrix.atoms)
    xyzarr = individual.matrix.converToXYZ()
    x, y, z = 0.0, 0.0, 0.0  # Note Unused
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
