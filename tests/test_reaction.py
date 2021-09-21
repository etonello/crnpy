#!/usr/bin/env python

"""Tests for Reaction class."""

from os import path
import unittest


from crnpy.crn import CRN, from_react_file
from crnpy.crncomplex import Complex
from crnpy.parsereaction import parse_reactions, parse_complex, parse_reaction, parse_expr
from crnpy.reaction import Reaction, translate
from .test_reduction import eqs_match

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


input_folder = path.join(path.dirname(__file__), 'data', 'reactions')

class TestReaction(unittest.TestCase):
    def test_parse_reaction(self):
        r, r_ = parse_reaction("     a + 2b <->(k1*a*b**2/(k2 +a)  ) 3   c       ")
        self.assertEqual(r.reactant, Complex({'a': 1, 'b': 2}))
        self.assertEqual(r.product, Complex({'c': 3}))
        self.assertEqual(r_.reactant, Complex({'c': 3}))
        self.assertEqual(r_.product, Complex({'a': 1, 'b': 2}))
        self.assertEqual(r.rate, parse_expr('k1*a*b**2/(k2+a)'))
        Reaction('r1', Complex(A = 1, B = 2), Complex(C = 1), parse_expr('k1*A*B**2'))

        r, r_ = parse_reactions(["a + 2b <->(k1/(k2+a)) 3c"])
        self.assertEqual(r.reactant, Complex({'a': 1, 'b': 2}))
        self.assertEqual(r.product, Complex({'c': 3}))
        self.assertEqual(r_.reactant, Complex({'c': 3}))
        self.assertEqual(r_.product, Complex({'a': 1, 'b': 2}))
        self.assertEqual(r.rate, parse_expr('k1*a*b**2/(k2+a)'))
        self.assertEqual(r_.kinetic_param, parse_expr('k_r0_rev'))

        r, r_ = parse_reactions(["a + 2b <->(k1/(k2+a)) 3c"], rate = True)
        self.assertEqual(r.rate, parse_expr('k1/(k2+a)'))
        self.assertEqual(r_.kinetic_param, parse_expr('k_r0_rev/c**3'))

        self.assertEqual(1, len(parse_reactions(['a -> b'])))
        self.assertRaises(ValueError, parse_reactions, ['a b'])
        self.assertRaises(ValueError, parse_reactions, 'a -> b')


    def test_rids(self):
        self.assertRaises(ValueError,
                          parse_reactions,
                          ["r1: a -> b", "c -> d", "e -> f", "r1: -> e"])
        self.assertEqual(['r1', 'r2', 'r0', 'r3'],
                [r.reactionid for r in parse_reactions(["r1: a -> b", "c -> d", "r0: e -> f", " -> e"])])


    def test_format(self):
        r = Reaction("r_123", Complex({"A": 1, "B": 2}), Complex({"B": 1, "C": 3}), parse_expr("k_123*A*B**2/(A+B+1)"))
        self.assertEqual((parse_expr("k_123/(A+B+1)") - parse_expr(r.format_kinetics())).cancel(), 0)
        self.assertEqual((parse_expr("k_123*A*B**2/(A+B+1)") - parse_expr(r.format_kinetics(rate = True))).cancel(), 0)
        r._kinetic_param = 0.00012345
        self.assertEqual(str(r), "r_123: A + 2B ->(1.234e-4) B + 3C")
        self.assertEqual(r.format(precision = 4), "r_123: A + 2B ->(1.2345e-4) B + 3C")
        self.assertEqual(r.format(precision = 1), "r_123: A + 2B ->(1.2e-4) B + 3C")
        self.assertEqual(r.format(True, 4), "r_123: A + 2B ->(1.2345e-4*A*B**2) B + 3C")
        self.assertEqual(r.format(True, precision = 1), "r_123: A + 2B ->(1.2e-4*A*B**2) B + 3C")


    def test_merge(self):
        crn = from_react_file(path.join(input_folder, "test_merge"))
        origspecies = crn.species
        origeqs = crn.equations()
        for r in crn.reactions: print(r)
        crn.merge_reactions()
        print
        for r in crn.reactions: print(r)
        eqs = crn.equations()
        couples = [(tuple(sorted(r.reactant.items())), tuple(sorted(r.product.items()))) for r in crn.reactions]
        self.assertEqual(len(couples), len(set(couples)))
        self.assertEqual(eqs_match(origeqs, origspecies, crn.removed_species, crn.equations(), crn.species), 0)


    def test_remove_react_prod(self):
        reaction = parse_reactions(["5 a + b + 2 c -> 3 a + 2 b + d"])[0]
        print(reaction)
        reaction.remove_react_prod()
        print(reaction)
        self.assertEqual(reaction.reactant, Complex({"a": 2, "c": 2}))
        self.assertEqual(reaction.product, Complex({"b": 1, "d": 1}))

        reaction = parse_reactions(["5 a + b + 2 c -> 3 a + 2 b + d"])[0]
        reaction.remove_react_prod('b')
        print(reaction)
        self.assertEqual(reaction.reactant, Complex({"a": 5, "c": 2}))
        self.assertEqual(reaction.product, Complex({"a":3, "b": 1, "d": 1}))


    def test_fix_ma(self):
        reaction = parse_reactions(["a ->(k1 * a * b**2 / (a + b**3)) c"])[0]
        print(reaction)
        print
        reaction._fix_ma(['a', 'b', 'c'])
        print(reaction)
        self.assertEqual(reaction.reactant, Complex({"a": 2, "b": 2}))
        self.assertEqual(reaction.product, Complex({"a": 1, "b": 2, "c": 1}))


    def test_fix_denom(self):
        reaction = parse_reactions(["a + 2 b + c + d ->(k1 / ((a + b**3)*c*b**2*d)) 2 b + 2 c"])[0]
        print(reaction)
        print
        reaction._fix_denom(['a', 'b', 'c'])
        print(reaction)
        self.assertEqual(reaction.reactant, Complex({"a": 1, "d": 1}))
        self.assertEqual(reaction.product, Complex({"c": 1}))


    def test_translate(self):
        r = Reaction('', parse_complex('x1 + x2'), parse_complex('2x1 + x3'), parse_expr('k*x1*x2'))
        r1 = Reaction('', parse_complex('4x1 + x2 + x4'), parse_complex('5x1 + x3 + x4'), parse_expr('k*x1*x2'))
        r2 = Reaction('', parse_complex('x2 + x4'), parse_complex('x1 + x3 + x4'), parse_expr('k*x1*x2'))
        self.assertEqual(r2, translate(r, Complex({'x1': -1, 'x4': 1})))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestReaction)
    unittest.TextTestRunner(verbosity=2).run(suite)
