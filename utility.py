from math import log
from random import random

class SSA:
    """Container for SSAs"""

    def __init__(self, model, max_t):
        """Initialize container with model and pseudorandom number generator"""
        self.model = model
        self.random = random
        self.max_t = max_t
    def direct(self):
        """Indefinite generator of direct-method trajectories"""
        while True:
            while not self.model.exit() and self.model["time"][-1] < self.max_t:
                # evaluate weights and partition
                weights = [
                    (rxn, sto, pro(self.model))
                    for (rxn, sto, pro) in self.model.reactions
                ]
                partition = sum(w[-1] for w in weights)

                # evaluate sojourn time (MC step 1)
                sojourn = log(1.0 / self.random()) / partition
                self.model["time"].append(self.model["time"][-1] + sojourn)

                # evaluate the reaction (MC step 2)
                partition = partition * self.random()
                while partition >= 0.0:
                    rxn, sto, pro = weights.pop(0)
                    partition -= pro
                for species, delta in sto.items():
                    self.model[species].append(self.model[species][-1] + delta)

                self.model.curate()
            yield self.model
            self.model.reset()

    def first_reaction(self):
        """Indefinite generator of 1st-reaction trajectories"""
        while True:
            while not self.model.exit() and self.model["time"][-1] < self.max_t:

                # evaluate next reaction times
                times = [
                    (
                        log(
                            1.0 / self.random()
                        ) / pro(self.model),
                        sto
                    )
                    for (rxn, sto, pro) in self.model.reactions
                ]
                times.sort()

                # evaluate reaction time
                self.model["time"].append(
                    self.model["time"][-1] + times[0][0]
                )

                # evaluate reaction
                for species, delta in times[0][1].items():
                    self.model[species].append(
                        self.model[species][-1] + delta
                    )

                self.model.curate()
            yield self.model
            self.model.reset()


class SSAModel(dict):
    """Container for SSA model"""

    def __init__(
            self, initial_conditions, propensities, stoichiometry
    ):
        """Initialize model"""
        super().__init__(**initial_conditions)
        self.reactions = list()
        self.excluded_reactions = list()
        for reaction,propensity in propensities.items():
            if propensity(self) == 0.0:
                self.excluded_reactions.append(
                    (
                        reaction,
                        stoichiometry[reaction],
                        propensity
                    )
                )
            else:
                self.reactions.append(
                    (
                        reaction,
                        stoichiometry[reaction],
                        propensity
                    )
                )

    def exit(self):
        """Return True to break out of trajectory"""

        # return True if no more reactions
        if len(self.reactions) == 0: return True

        # return False if there are more reactions
        else: return False

    def curate(self):
        """Validate and invalidate model reactions"""

        # evaluate possible reactions
        reactions = []
        while len(self.reactions) > 0:
            reaction = self.reactions.pop()
            if reaction[2](self) == 0:
                self.excluded_reactions.append(reaction)
            else:
                reactions.append(reaction)
        self.reactions = reactions

        # evaluate impossible reactions
        excluded_reactions = []
        while len(self.excluded_reactions) > 0:
            reaction = self.excluded_reactions.pop()
            if reaction[2](self) > 0:
                self.reactions.append(reaction)
            else:
                excluded_reactions.append(reaction)
        self.excluded_reactions = excluded_reactions

    def reset(self):
        """Clear the trajectory"""

        # reset species to initial conditions
        for key in self: del self[key][1:]

        # reset reactions per initial conditions
        self.curate()


