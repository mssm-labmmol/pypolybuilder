
import numpy as np
from deap import base
from deap import creator
from deap import tools
import math


class PSO:
    def __init__(self, evaluateFuntion, indSize, weight=(-1.0,), popSize=100, rangeGen=[0, 1], smin = 0.25, smax = 0.75, phi1=2.0, phi2=2.0):

        self.weight = weight
        self.indSize = indSize
        self.popSize = int(popSize)
        self.rangeGen = rangeGen
        self.evaluateFuntion = evaluateFuntion
        self.smin = smin
        self.smax = smax
        self.phi1 = phi1
        self.phi2 = phi2


    def run(self, NGENS=100):

        creator.create("FitnessMax", base.Fitness, weights=self.weight)
        creator.create("Particle", np.ndarray, fitness=creator.FitnessMax, speed=list,smin=None, smax=None, best=None)

        def generate(size, pmin, pmax, smin, smax):
            part = creator.Particle(np.random.uniform(pmin, pmax, size))
            part.speed = np.random.uniform(smin, smax, size)
            part.smin = smin
            part.smax = smax
            return part

        def updateParticle(part, best, phi1, phi2):
            u1 = np.random.uniform(0, phi1, len(part))
            u2 = np.random.uniform(0, phi2, len(part))
            v_u1 = u1 * (part.best - part)
            v_u2 = u2 * (best - part)
            part.speed += v_u1 + v_u2
            for i, speed in enumerate(part.speed):
                if abs(speed) < part.smin:
                    part.speed[i] = math.copysign(part.smin, speed)
                elif abs(speed) > part.smax:
                    part.speed[i] = math.copysign(part.smax, speed)
            part += part.speed

        toolbox = base.Toolbox()
        toolbox.register("particle", generate, size=self.indSize, pmin=self.rangeGen[0], pmax=self.rangeGen[1], smin=self.smin, smax=self.smax)
        toolbox.register("population", tools.initRepeat, list, toolbox.particle)
        toolbox.register("update", updateParticle, phi1=self.phi1, phi2=self.phi2)
        toolbox.register("evaluate", self.evaluateFuntion)



        pop = toolbox.population(n=self.popSize);
        #hof = tools.HallOfFame(1)
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)


        logbook = tools.Logbook()
        logbook.header = "gen", "evals", "std", "min", "avg", "max"

        # Evaluate the individuals
        fitnesses = toolbox.map(toolbox.evaluate, pop)
        for ind, fit in zip(pop, fitnesses):
            ind.fitness.values = fit

        record = stats.compile(pop)
        logbook.record(gen=0, evals=len(pop), **record)
        #print(logbook.stream)

        best=None
        #hof.update(pop)
        for g in range(0, NGENS):
            for part in pop:
                part.fitness.values = toolbox.evaluate(part)
                if part.best is None or part.best.fitness < part.fitness:
                    part.best = creator.Particle(part)
                    part.best.fitness.values = part.fitness.values
                if best is None or best.fitness < part.fitness:
                    best = creator.Particle(part)
                    best.fitness.values = part.fitness.values
            for part in pop:
                toolbox.update(part, best)

            logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
            #print(logbook.stream)
        return best, best.fitness.values[0]