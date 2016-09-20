#!/usr/bin/env python

"""Tests for reduction methods of CRN class.
References
==========
[1] Cornish-Bowden, Athel. "Fundamentals of enzyme kinetics". Elsevier, 2014.
[2] Ingalls, Brian. "Mathematical Modelling in Systems Biology: An Introduction.", 2013.
[3] SBMLsqueezer documentation, http://www.cogsys.cs.uni-tuebingen.de/software/SBMLsqueezer/doc/KineticLaws2.pdf.
[4] Segel, Irwin H. "Enzyme kinetics". Vol. 957. Wiley, New York, 1975.
"""

from os import listdir, path
import sympy as sp
import sys
import unittest
import warnings


from crnpy.conslaw import ConsLaw
from crnpy.crn import CRN, from_sbml, from_react_file, from_react_strings
from crnpy.crncomplex import to_complex, sympify
from crnpy.parsereaction import parse_reactions

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


input_reactions = path.join(path.dirname(__file__), 'data', 'reactions')
input_sbml = path.join(path.dirname(__file__), 'data', 'sbml')


class TestReduction(unittest.TestCase):
    def test_enzyme(self):
        """One-substrate enzyme kinetics."""

        filename = path.join(input_sbml, "enzyme.xml")
        crn = from_sbml(filename)
        crn.save_sbml(path.join(input_sbml, "enzyme_original.xml"))
        crn = from_sbml(path.join(input_sbml, "enzyme_original.xml"))

        rate = sympify("_comp*E*vcat_kcat*veq_kon*S/(vcat_kcat + veq_koff)")

        crn.qss('ES')
        self.assertEqual((crn.rates[0] - rate).simplify(), 0)
        crn.save_sbml(path.join(input_sbml, "enzyme_simplified_qss_with_e.xml"))

        crn.remove_constant('E')
        self.assertEqual((crn.rates[0] - rate).simplify(), 0)
        crn.save_sbml(path.join(input_sbml, "enzyme_simplified_qss.xml"))

        # Michaelis-Menten
        crn = from_sbml(filename)
        enzyme_cons_law = ConsLaw('E + ES', 'Et')
        crn.qss(cons_law = ('E', enzyme_cons_law))
        rate = sympify("_comp*Et*vcat_kcat*veq_kon*S/(vcat_kcat + veq_koff + veq_kon*S)")
        self.assertEqual((crn.rates[0] - rate).simplify(), 0)
        crn.save_sbml(path.join(input_sbml, "enzyme_simplified_MM.xml"))
        crn.save_reaction_file(path.join(input_reactions, "enzyme_simplified_MM"))

        # Reading Michaelis-Menten model
        crn = from_sbml(path.join(input_sbml, "enzyme_simplified_MM.xml"))


    def test_enzyme_rapid_eq(self):
        """One-substrate enzyme kinetics: rapid equilibrium."""
        # Rapid equilibrium
        filename = path.join(input_sbml, "enzyme.xml")
        crn = from_sbml(filename)
        enzyme_cons_law = ConsLaw('E + ES', 'Et')
        crn.rapid_eq('ES', 'S + E', cons_law = ('E', enzyme_cons_law))
        self.assertEqual((crn.kinetic_params[0] - sympify("Et*vcat_kcat*_comp/(S + veq_koff/veq_kon)")).simplify(), 0)


    def test_enzyme_reversible(self):
        """One-substrate enzyme reversible kinetics."""
        crn = from_react_file(path.join(input_reactions, "enzyme_reversible"))

        crn.qss(cons_law = ('E', ConsLaw('E + C', 'Et')))
        forward = sympify("Vf*S/Ks/(1 + S/Ks + P/Kp)").subs(sympify("Vf"), sympify("Et*k2"))\
                                                         .subs(sympify("Ks"), sympify("(k_1+k2)/k1"))\
                                                         .subs(sympify("Kp"), sympify("(k_1+k2)/k_2")).simplify()
        backward = sympify("Vr*P/Kp/(1 + S/Ks + P/Kp)").subs(sympify("Vr"), sympify("Et*k_1"))\
                                                          .subs(sympify("Ks"), sympify("(k_1+k2)/k1"))\
                                                          .subs(sympify("Kp"), sympify("(k_1+k2)/k_2")).simplify()
        self.assertEqual(set(crn.rates), set([forward, backward]))


    def test_multi_product(self):
        """Test for the case of intermediate produced
        with multiple stoichiometries."""
        crn = from_react_file(path.join(input_reactions, "multi_product"))
        origspecies = crn.species
        origeqs = crn.equations()
        self.assertEqual(eqs_match(origeqs, origspecies, crn.removed_species, crn.equations(), crn.species), 0)
        crn._qss_generalised('y', no_rates = True)
        reacts = parse_reactions(["r0_3r1: a + 2d ->(k_r0_3r1) 3b + c",
                                  "r0_3r4: a ->(k_r0_3r4) 3a + 6b + c + d",
                                  "r2_r1: d + e + f ->(k_r2_r1) b",
                                  "r2_r4: e + f ->(k_r2_r4) a + 2b",
                                  "r3_2r1: 2d + 2e ->(k_r3_2r1) 2b + h",
                                  "r3_2r4: 2e ->(k_r3_2r4) 2a + 4b + h"])
        self.assertTrue(all(r in reacts for r in crn.reactions) and all(r in crn.reactions for r in reacts))


    def test_examples(self):
        """Examples."""
        folder = path.join(input_reactions, "examples/")
        files = [f for f in listdir(folder) if path.isfile(path.join(folder, f))]
        solFolder = path.join(input_reactions, "examples/sols")
        notVerified = []

        for reactionFile in files:
            print("######################################")
            print("Example: {}".format(reactionFile))
            crn = from_react_file(path.join(folder, reactionFile))
            origspecies = crn.species
            origeqs = crn.equations()
            crn.qss('i')

            # Check that equations are as expected
            self.assertEqual(eqs_match(origeqs, origspecies, crn.removed_species, crn.equations(), crn.species), 0)

            # Check that reactions are as in solution file
            reactions = map(lambda x: x.split(': ')[1], map(str, crn.reactions))
            solFile = path.join(solFolder, reactionFile + "_sol")
            if path.isfile(solFile):
                with open(solFile) as f:
                    reacts = map(lambda x: x.split(': ')[1], f.read().splitlines())
                self.assertEqual(set(reactions), set(reacts))
            else: notVerified.append(reactionFile)

        if len(notVerified) > 0: print("Solutions not checked: {}".format(notVerified))


    def test_qss1(self):
        """QSS test 1 (Ingalls, section 2.2.1)."""
        crn = from_react_file(path.join(input_reactions, "basic1"))
        crn.qss('b')
        inda = crn.complexes.index(to_complex('a'))
        self.assertEqual((crn.laplacian[inda,inda]-sympify("k1*k2/(k2 + k_1)")).simplify(), 0)


    def test_qss2(self):
        """QSS test 2 (Ingalls, section 2.2.1)."""
        crn = from_react_file(path.join(input_reactions, "basic2"))
        crn.qss('a')
        rates = map(sympify, ["k0", "k2 * b"])
        self.assertEqual(set(rates), set(crn.rates))


    def test_qss3(self):
        """QSS test 3 (Ingalls, section 2.2.1)."""
        crn = from_react_file(path.join(input_reactions, "basic3"))
        crn.qss('a')
        indb = crn.complexes.index(to_complex('b'))
        self.assertEqual((crn.laplacian[indb,indb]-sympify("(k0*k_1/(k0 + k1)+k2)")).simplify(), 0)


    def test_rapid_eq1(self):
        """Rapid equilibrium with pooling test 1 (Ingalls, section 2.2.1)."""
        crn = from_react_file(path.join(input_reactions, "basic1"))
        crn.rapid_eq_with_pool('b', 'a', pool_name = "c")
        self.assertEqual(sp.simplify(crn.laplacian[0,0] - sympify("k2 * k1 / (k_1 + k1)")), 0)


    def test_rapid_eq2(self):
        """Rapid equilibrium with pooling test 2 (Ingalls, section 2.2.1)."""
        crn = from_react_file(path.join(input_reactions, "basic2"))
        crn.rapid_eq_with_pool('b', 'a', pool_name = "c")
        rates = map(sympify, ["k0", "k2 * k1 * c / (k_1 + k1)"])
        self.assertEqual(set(rates), set(crn.rates))


    def test_rapid_eq3(self):
        """Rapid equilibrium with pooling test 3 (Ingalls, exercise 2.2.1)."""
        crn = from_react_file(path.join(input_reactions, "basic3"))
        crn.rapid_eq_with_pool('b', 'a', pool_name = "c")
        self.assertEqual(sp.simplify(crn.laplacian[0,0] - sympify("(k0 * k_1 + k2 * k1)/(k_1 + k1)")), 0)


    def test_rapid_eq4(self):
        reacts = ["2 a + b ->(kf) y",
                  "y ->(kr) 2 a + b",
                  "y ->(k) d",
                  "y + d ->(k1) a",
                  "a ->(k2) y + d"]
        crn = from_react_strings(reacts)
        crn.remove(rapid_eq = [('y', '2 * a + b')])
        sol = parse_reactions(["2a + b ->(k*kf/kr) d", "2a + b + d (k2)<->(k1*kf/kr) a"])
        self.assertTrue(all(r in crn.reactions for r in sol))
        self.assertTrue(all(r in sol for r in crn.reactions))

        reacts = ["a + 2 y ->(k1) b",
                  "3a + c ->(k2) y",
                  "y ->(k_2) 3a + c",
                  "y ->(k3) d"]
        crn = from_react_strings(reacts)
        crn.remove(rapid_eq = [('y', '3*a + c')])
        sol = parse_reactions(["7a + 2c ->(k1*k2^2/k_2^2) b", "3a + c ->(k3*k2/k_2) d"])
        self.assertTrue(all(r in crn.reactions for r in sol))
        self.assertTrue(all(r in sol for r in crn.reactions))


    def test_qss_errors(self):
        crn = from_react_strings(['a -> b', 'b + d ->(k*b) c + d', 'd -> '])
        self.assertRaises(ValueError, crn.qss, 'a') # a not intermediate
        self.assertRaises(ValueError, crn.qss, 'b') # nonlinear kinetics
        self.assertRaises(ValueError, crn.qss, 'c') # c not intermediate
        self.assertRaises(ValueError, crn.qss, 'd') # d not intermediate
        self.assertRaises(ValueError, crn.qss, 'e') # e not in valid species


    def test_re_errors(self):
        crn = from_react_strings(['a (k_0)<->(k0) b + c', 'b (k_1)<->(k1*c) c', 'd ->(k2) '])
        self.assertRaises(TypeError, crn.rapid_eq, 'a') # invalid input
        self.assertRaises(ValueError, crn.rapid_eq, 'a', 'f') # f not valid complex
        self.assertRaises(ValueError, crn.rapid_eq, 'a + c', 'b') # a + c not valid species
        self.assertRaises(ValueError, crn.rapid_eq, 'c', 'b') # unable to remove c
        crn = from_react_strings(['a (k_0)<->(k0) b + c', 'b (k_1)<->(k1*b) c', 'd ->(k2) '])
        self.assertRaises(ValueError, crn.rapid_eq, 'b', 'c') # unable to remove b


    def test_qss_minimal(self):
        reacts = ["a -> y + b", "y + b -> z", "z -> y", "y -> d"]

        crn = from_react_strings(reacts)
        # Removing z first
        crn.qss('z')
        crn.qss('y')
        reactions = crn.reactions

        crn = from_react_strings(reacts)
        # Removing y first
        crn.qss('y')
        crn.qss('z')
        crn = from_react_strings(reacts)
        # Removing using minimal
        crn.remove(qss = ['y', 'z'], minimal = True)

        self.assertEqual(len(reactions), len(crn.reactions))
        self.assertTrue(all(crn.reactions[i] in reactions for i in range(crn.n_reactions)))


    def test_qss_minimal2(self):
        reacts = ["r1: x -> y",
                  "r1_rev: y -> x",
                  "r2: e + x -> w",
                  "r2rev: w -> e + x",
                  "r3: f + y -> z",
                  "r3rev: z -> f + y",
                  "r4: w -> z",
                  "r5: d + x -> c",
                  "r5rev: c -> d + x",
                  "r6: b + y -> a",
                  "r6rev: a -> b + y"]

        crn = from_react_strings(reacts)
        crn.set_params(dict((k, 1) for k in crn.kinetic_params))
        # Qss on y, z, w, x
        crn.remove(qss = ['y', 'z', 'w', 'x'], minimal = True)
        self.assertTrue(len(crn.reactions) == 4)
        for r in crn.reactions: print(r)
        print

        crn = from_react_strings(reacts)
        crn.set_params(dict((k, 1) for k in crn.kinetic_params))
        # Qss on y, z, x, w
        crn.remove(qss = ['y', 'z', 'x', 'w'], minimal = True)
        self.assertTrue(len(crn.reactions) == 4)
        for r in crn.reactions: print(r)


def eqs_match(origeqs, origspecies, removed, neweqs, newspecies, debug = False):
    """Check that the original equations match the new equations after
    replacing the species with the expression in removed."""
    origeqsdict = dict((origspecies[s], origeqs[s]) for s in range(len(origspecies)))
    replaceeqs = {}
    for s in origeqsdict:
        e = origeqsdict[s]
        for i in range(len(removed)): e = e.subs(removed[i][0], removed[i][1])
        e.simplify()
        replaceeqs[s] = e
    neweqsdict = dict((newspecies[s], neweqs[s]) for s in range(len(newspecies)))

    if debug:
        for s in newspecies: print("{}: {}".format(s, replaceeqs[s] - neweqsdict[s]).simplify())

    return sum((replaceeqs[s] - neweqsdict[s]).simplify() != 0 for s in newspecies)


if __name__ == "__main__":
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestReduction)
    unittest.TextTestRunner(verbosity=2).run(suite)
