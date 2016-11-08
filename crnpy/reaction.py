#!/usr/bin/env python

"""Reaction class and functions."""

from collections import defaultdict
from functools import reduce
from libsbml import formulaToL3String, parseL3Formula
import sympy as sp
import copy

from .crncomplex import Complex

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


class Reaction(object):
    """A reaction is defined by a string reactionid,
    a reactant complex, a product complex,
    and a string or sympy expression representing the rate.

    :Example:

    >>> from crnpy.reaction import Reaction
    >>> from crnpy.crncomplex import Complex
    >>> r = Reaction("r1", Complex(A = 1, B = 2), Complex(C = 1), "k1*A*B**2")
    >>> r
    r1: A + 2B ->(k1) C
    >>> r.reactionid, r.reactant, r.product, r.rate, r.kinetic_param
    ('r1', A + 2B, C, A*B**2*k1, k1)

    Attributes: reactionid, reactant, product, rate, kinetic_param.
    """
    def __init__(self, reactionid, reactant, product, rate):
        self._reactionid = reactionid
        self._reactant = reactant
        self._product = product
        self._rate = rate

    def __copy__(self):
        return Reaction(self.reactionid, self.reactant, self.product, self.rate)
    def __deepcopy__(self, memo):
        return Reaction(copy.deepcopy(self.reactionid, self.reactant, self.product, self.rate, memo))

    @property
    def reactionid(self):
        """String id of the reaction.

        :type: string.
        """
        return self._reactionid

    @property
    def reactant(self):
        """Reactant complex of the reaction.

        :type: Complex.
        """
        return self._reactant

    @property
    def product(self):
        """Product complex of the reaction.

        :type: Complex.
        """
        return self._product

    @property
    def rate(self):
        """Rate of the reaction.

        :type: sympy expression.
        """
        return self._rate

    @property
    def _rate(self):
        return self.__rate
    @_rate.setter
    def _rate(self, value):
        self.__rate = value
        if value is not None: self.__kinetic_param = self.rate / self.reactant.ma()

    @property
    def kinetic_param(self):
        """Rate of the reaction divided by the mass action monomial
        associated to the reactant.

        :type: sympy expression.
        """
        return self._kinetic_param

    @property
    def _kinetic_param(self):
        return self.__kinetic_param
    @_kinetic_param.setter
    def _kinetic_param(self, value):
        self._rate = (value * self.reactant.ma()).cancel()

    def __str__(self):
        return self.format()

    def  __repr__(self):
        return self.__str__()

    def __eq__(self, reaction):
        return self.reactant == reaction.reactant and \
               self.product == reaction.product and \
               (((self.rate - reaction.rate) == 0) or
                ((self.rate - reaction.rate).cancel() == 0))

    def format(self, rate = False, precision = 3):
        """Return a string of the form
        reactant complex ->(k) product complex
        where k is the generalised kinetic parameter of the reaction if rate = False,
        otherwise the rate of the reaction.

        :rtype: string.
        """
        return "{}: {} ->{} {}".format(self.reactionid,
                                       self.reactant,
                                       "(" + self.format_kinetics(rate, precision) + ")" if self.rate else "",
                                       self.product)

    def format_kinetics(self, rate = False, precision = 3):
        """Convert the kinetic parameter or rate to string.
        If rate = True, return a string representing the rate,
        otherwise one representing the kinetic parameter.
        If the kinetic parameter is a float, use exponent notation,
        with a number of digits equal to precision (3 is the default).
        If the reaction has no defined rate, return None."""
        k = None
        if self.rate:
            if isinstance(self.kinetic_param, sp.Float):
                k = "{:.{}e}".format(self.kinetic_param, precision)
                if rate:
                    k = k + "*" + str(self.reactant.ma())
            else:
                if rate:
                    k = str(self.rate)
                else:
                    k = str(self.kinetic_param)
        return k

    def latex(self, rate = False):
        """Return the latex code for the reaction.
        By default the kinetic parameter of the reaction is included.
        To use the rate instead, use rate = True.

        :Example:

        >>> from crnpy.reaction import Reaction
        >>> from crnpy.crncomplex import Complex
        >>> r = Reaction("r1", Complex(A = 1, B = 2), Complex(C = 1), "k1*A*B**2")
        >>> print(r.latex())
        r_{1}: A + 2 B \\xrightarrow{k_{1}} C
        >>> print(r.latex(True))
        r_{1}: A + 2 B \\xrightarrow{A B^{2} k_{1}} C

        :rtype: string.
        """
        return "{}: {} {} {}".format(sp.latex(self.reactionid),
                                     sp.latex(self.reactant.symp()),
                                     str("\\xrightarrow{" + sp.latex(self.rate if rate else self.kinetic_param) + "}") if self.rate else str("\\rightarrow"),
                                     sp.latex(self.product.symp()))

    def remove_react_prod(self, species = None):
        """Remove common species between reactant and product.

        If a species is specified, only that species is removed.

        :Example:

        >>> reacts = parse_reactions(["A + 2B -> A + B", "A + B + C -> B + C + D"])
        >>> reacts[0].remove_react_prod()
        >>> reacts[1].remove_react_prod('C')
        >>> reacts
        [r0: B ->(A*B*k_r0) , r1: A + B ->(C*k_r1) B + D]

        """
        reactant = Complex(self.reactant)
        product = Complex(self.product)
        if species == None:
            common = reactant & product
            self._reactant = Complex(reactant - common)
            self._product = Complex(product - common)
        else:
            if species in reactant and species in product:
                if reactant[species] > product[species]:
                    self.reactant[species] = reactant[species] - product[species]
                    del self.product[species]
                if product[species] > reactant[species]:
                    self.product[species] = product[species] - reactant[species]
                    del self.reactant[species]
                if reactant[species] == product[species]:
                    del self.reactant[species]
                    del self.product[species]
        # Adjust kinetic parameter
        self.__kinetic_param = (self.rate / self.reactant.ma()).cancel()


    def _fix_ma(self, species = None):
        """Check the numerator of the reaction rate, and adds species
        to reactant and product if they divide the numerator but their
        stoichiometry does not match the degree in the rate."""
        remainder = self.kinetic_param.as_numer_denom()[0].cancel()

        if remainder.func.__name__ == 'Mul':
            mulargs = list(remainder.args) + [i.args[0] for i in remainder.args if i.func.__name__ == 'Mul'] \
                                           + [i.args[0] for i in remainder.args if i.func.__name__ == 'Pow']
            while any(sp.Symbol(s) in mulargs for s in species):
                for s in species:
                    if sp.Symbol(s) in mulargs:
                        if s in self.reactant: self.reactant[s] = self.reactant[s] + 1
                        else: self.reactant[s] = 1
                        if s in self.product: self.product[s] = self.product[s] + 1
                        else: self.product[s] = 1
                        remainder = (remainder / sp.Symbol(s)).factor()
                        if remainder.func.__name__ == 'Mul':
                            mulargs = list(remainder.args) + [i.args[0] for i in remainder.args if i.func.__name__ == 'Mul'] \
                                                           + [i.args[0] for i in remainder.args if i.func.__name__ == 'Pow']
                        else: mulargs = []
            # update the kinetic parameter
            self.__kinetic_param = (self.rate / self.reactant.ma()).cancel()


    def _fix_denom(self, species):
        """Remove species that are involved in both reactant and product,
        if their concentration divides the denominator of the rate."""
        remainder = self.kinetic_param.as_numer_denom()[1].cancel()

        #if remainder.func.__name__ == 'Mul':
        if remainder != 1:
            mulargs = [remainder] + list(remainder.args) + [i.args[0] for i in remainder.args if i.func.__name__ == 'Mul'] \
                                                         + [i.args[0] for i in remainder.args if i.func.__name__ == 'Pow']
            while any(sp.Symbol(s) in mulargs and s in self.reactant and s in self.product for s in species):
                for s in species:
                    if sp.Symbol(s) in mulargs and s in self.reactant and s in self.product:
                        if self.reactant[s] == 1: del self.reactant[s]
                        else: self.reactant[s] = self.reactant[s] - 1
                        if self.product[s] == 1: del self.product[s]
                        else: self.product[s] = self.product[s] - 1
                        remainder = (remainder / sp.Symbol(s)).factor()
                        if remainder.func.__name__ == 'Mul':
                            mulargs = list(remainder.args) + [i.args[0] for i in remainder.args if i.func.__name__ == 'Mul'] \
                                                           + [i.args[0] for i in remainder.args if i.func.__name__ == 'Pow']
                        else:
                            if str(remainder) in species: mulargs = [remainder]
                            else: mulargs = []
        # update the kinetic parameter
        self._kinetic_param = self.rate / self.reactant.ma()


def _split_reaction(reaction):
    """Split a reaction into separate reactions, one
    for each addend in rate (compare to split_reactions_monom).

    :rtype: list of Reactions.
    """
    ratenumer, ratedenom = reaction.rate.as_numer_denom()
    ratenumer = ratenumer.expand()
    if ratenumer.func.__name__ == 'Add':
        reactions = []
        rateadds = list(ratenumer.args)

        for ra in range(len(rateadds)):
            reactions.append(Reaction(reaction.reactionid + "_" + str(ra), \
                                      reaction.reactant, \
                                      reaction.product, \
                                      rateadds[ra] / ratedenom))
        return reactions
    else: return [reaction]


def _split_reaction_monom(reaction, species):
    """Split a reaction into separate reactions, one
    for each monomial in rate (compare to split_reactions).

    :rtype: list of Reactions.
    """
    ratenumer, ratedenom = reaction.rate.cancel().as_numer_denom()
    ratenumer = ratenumer.expand()
    species = map(sp.Symbol, species)
    ratendict = sp.Poly(ratenumer, *species).as_dict()
    if len(ratendict) > 1:
        reactions = []

        i = 0
        for degrees in ratendict:
            i = i + 1
            ratenpart = sp.Mul(*[species[r]**degrees[r] for r in range(len(species))]) * ratendict[degrees]
            reactions.append(Reaction(reaction.reactionid + "_" + str(i), \
                                      reaction.reactant, \
                                      reaction.product, \
                                      ratenpart / ratedenom))
        return reactions
    return [reaction]


def merge_reactions(reactions):
    """Merge reactions with same reactants and products.

    Take a list of reactions in input and return a list of reactions,
    with a maximum of one reaction for each pair of reactant and product.
    The rates of reactions with the same reactant and product are summed, and
    their reaction ids are concatenated.

    :Example:

    >>> from crnpy.parsereaction import parse_reactions
    >>> from crnpy.reaction import merge_reactions
    >>> reacts = parse_reactions(["a -> b", "c <-> d + e", "d + e -> c", "a -> b"])
    >>> merge_reactions(reacts)
    [r0r3: a ->(k_r0 + k_r3) b, r1: c ->(k_r1) d + e, r1_revr2: d + e ->(k_r1_rev + k_r2) c]

    :rtype: list of Reactions.
    """
    react = defaultdict(list)
    newreactions = []
    for reaction in reactions:
        react[(tuple(sorted(reaction.reactant.items())), tuple(sorted(reaction.product.items())))].append(reaction)
    for c in react:
        if react[c][0].reactant != react[c][0].product:
            newreactions.append(Reaction(''.join([reaction.reactionid for reaction in react[c]]), \
                                         react[c][0].reactant, \
                                         react[c][0].product, \
                                         sum([reaction.rate for reaction in react[c]]).factor()))
    return sorted(newreactions, key = lambda r: r.reactionid)


def _same_denom(reactions):
    """Change the rates so that they all have the same denominator.

    :rtype: list of Reactions.
    """
    numers, denoms = zip(*[reaction.rate.as_numer_denom() for reaction in reactions])
    commondenom = reduce(sp.lcm, denoms)
    newreactions = []
    for r in range(len(reactions)):
        reaction = reactions[r]
        if (denoms[r] - commondenom).cancel() != 0:
            diff = (commondenom / denoms[r]).cancel().expand()
            if diff.func.__name__ == 'Add':
                rateadds = list(diff.args)
                for ra in range(len(rateadds)):
                    newreactions.append(Reaction(reaction.reactionid + "_" + str(ra), \
                                                 reaction.reactant, \
                                                 reaction.product, \
                                                 rateadds[ra] * numers[r] / commondenom))
            else:
                newreactions.append(Reaction(reaction.reactionid, \
                                             reaction.reactant, \
                                             reaction.product, \
                                             diff * numers[r] / commondenom))
        else: newreactions.append(reactions[r])
    return newreactions


def translate(reaction, c):
    """Translate the reaction by c.
    Return the reaction (r + Id_c), where c has been added to both reactant and product.
    The rate is unchanged.

    :param reaction: reaction to translate.
    :type reaction: Reaction.
    :param c: complex to add to reactant and product.
    :type c: string.
    :rtype: Reaction.
    """
    rid = reaction.reactionid + "_" + str(c).replace(" ", "").replace("+", "_").replace("-", "m")
    return Reaction(rid, Complex(reaction.reactant + c), Complex(reaction.product + c), reaction.rate)


def reaction_path(reactions):
    """Translate the reactions so that they can be composed,
    in the given order.

    :rtype: list of Reactions.
    """
    additions = [reduce(lambda a, b: a + b,
                        [reactions[i].product for i in range(h)] +
                        [reactions[i].reactant for i in range(h+1, len(reactions))])
                 for h in range(len(reactions))]
    c = Complex(reduce(lambda a, b: a & b, additions))
    additions = [Complex(additions[h]-c) for h in range(len(additions))]
    return [translate(reactions[h], additions[h]) for h in range(len(reactions))], additions
