#!/usr/bin/env python

"""Tests for CRN class."""

import unittest

from filecmp import cmp
import libsbml
import sympy as sp
from os import path

from context import crnpy

from crnpy.conslaw import ConsLaw
from crnpy.createmodel import model_from_reacts
from crnpy.crn import CRN, from_sbml, from_react_file, from_reacts, from_react_strings
from crnpy.crncomplex import Complex, to_complex, sympify
from crnpy.parsereaction import parse_reaction_file, parse_reactions
from crnpy.reaction import Reaction

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


input_reactions = path.join(path.dirname(__file__), 'data', 'reactions')
input_sbml = path.join(path.dirname(__file__), 'data', 'sbml')


class TestCrn(unittest.TestCase):
    def test_empty_crn(self):
        for crn in [CRN(), from_reacts([]), from_react_strings([])]:
            self.assertEqual(0, crn.n_species)
            self.assertEqual((), crn.species)
            self.assertEqual(0, crn.n_complexes)
            self.assertEqual((), crn.complexes)
            self.assertEqual(0, crn.n_reactions)
            self.assertEqual((), crn.reactions)
            self.assertEqual((), crn.cons_laws)
            self.assertEqual((), crn.removed_species)
            self.assertTrue(crn.is_ma)
            self.assertEqual(sp.Matrix(), crn.complex_matrix)
            self.assertEqual(sp.Matrix(), crn.incidence_matrix)
            self.assertEqual(sp.Matrix(), crn.stoich_matrix)
            self.assertEqual(sp.Matrix(), crn.laplacian)
            self.assertEqual(sp.Matrix(), crn.influence_matrix())


    def test_empty_replace(self):
        crn = CRN()
        crn.reactions = parse_reactions(['a -> b', 'c + 2d <-> e', 'e -> a'])
        self.assertEqual(5, crn.n_species)
        self.assertEqual(('a', 'b', 'c', 'd', 'e'), crn.species)
        self.assertEqual(4, crn.n_complexes)
        self.assertEqual(['a', 'b', 'c + 2d', 'e'], list(map(str, crn.complexes)))
        self.assertEqual(4, crn.n_reactions)
        self.assertEqual(None, crn.model)
        crn.update_model()
        self.assertNotEqual(None, crn.model)
        self.assertEqual(['a', 'b', 'c', 'd', 'e'], [crn.model.getSpecies(s).getId() for s in range(crn.model.getNumSpecies())])
        self.assertEqual(4, crn.model.getNumReactions())
        Y = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 2, 0], [0, 0, 0, 1]])
        self.assertEqual(Y, crn.complex_matrix)
        Ia = sp.Matrix([[-1, 0, 0, 1], [1, 0, 0, 0], [0, -1, 1, 0], [0, 1, -1, -1]])
        self.assertEqual(Ia, crn.incidence_matrix)
        S = sp.Matrix([[-1, 0, 0, 1], [1, 0, 0, 0], [0, -1, 1, 0], [0, -2, 2, 0], [0, 1, -1, -1]])
        self.assertEqual(S, crn.stoich_matrix)
        k_r0, k_r1, k_r1_rev, k_r2 = crn.kinetic_params
        L = sp.Matrix([[k_r0, 0, 0, -k_r2], [-k_r0, 0, 0, 0], [0, 0, k_r1, -k_r1_rev], [0, 0, -k_r1, k_r1_rev + k_r2]])
        self.assertEqual(L, crn.laplacian)
        g_1_1, g_3_2, g_4_2, g_5_3, g_5_4 = sp.symbols('g_1_1 g_3_2 g_4_2 g_5_3 g_5_4')
        im = sp.Matrix([[g_1_1, 0, 0, 0], [0, 0, 0, 0], [0, g_3_2, 0, 0], [0, g_4_2, 0, 0], [0, 0, g_5_3, g_5_4]])
        self.assertEqual(im, crn.influence_matrix())


    def test_crn_from_reacts(self):
        reactions = [
                ("r1", {"A": 1, "B": 1}, {"C": 1}, "k1 * A * B"),
                ("r2", {"A": 1, "B": 1}, {"D": 1}, "k2 * A * B"),
                ("r3", {"C": 1}, {"E": 1, "F": 1}, "k3 * C"),
                ("r4", {"C": 1, "D": 1}, {"A": 1, "D": 1}, "k4 * C * D")]
        filename = path.join(input_sbml, "test_model1.xml")

        model, document, _ = model_from_reacts(list(map(lambda x: Reaction(x[0], \
                                                                  Complex(x[1]), \
                                                                  Complex(x[2]), \
                                                                  rate = sp.sympify(x[3])), reactions)))
        success = libsbml.writeSBMLToFile(document, filename)
        self.assertTrue(success)

        crn = from_sbml(filename)
        crn.inspect(True)


    def test_crn_from_react_file(self):
        reactions = parse_reaction_file(path.join(input_reactions, "allosteric_activation"))
        filename = path.join(input_sbml, "allosteric_activation.xml")

        model, document, _ = model_from_reacts(reactions)
        success = libsbml.writeSBMLToFile(document, filename)
        self.assertTrue(success)

        crn = from_sbml(filename)
        crn.inspect()


    def test_crn_from_sbml(self):
        filename = path.join(input_sbml, "BIOMD0000000001.xml")
        output = path.join(input_sbml, "out_BIOMD0000000001.xml")

        crn = from_sbml(filename)
        success = libsbml.writeSBMLToFile(crn.document, output)
        self.assertTrue(success)


    def test_reaction_setter(self):
        net = CRN()
        net.reactions = parse_reactions(["A -> 2B", "B <-> C + D"])
        self.assertEqual(("A", "B", "C", "D"), net.species)
        self.assertEqual(4, net.n_species)
        self.assertEqual(4, net.n_complexes)
        self.assertEqual(3, net.n_reactions)
        self.assertEqual(("r0", "r1", "r1_rev"), net.reactionids)
        self.assertEqual(sorted(["A", "2B", "B", "C + D"]), sorted(list(map(str, net.complexes))))
        self.assertEqual(("r0: A ->(k_r0) 2B", "r1: B ->(k_r1) C + D", "r1_rev: C + D ->(k_r1_rev) B"), tuple(map(str, net.reactions)))
        self.assertEqual(sp.Matrix([sympify("k_r0*A"), sympify("k_r1*B"), sympify("k_r1_rev*C*D")]), net.rates)
        S = sp.Matrix([[-1,  0,  0], [ 2, -1,  1], [ 0,  1, -1], [ 0,  1, -1]])
        self.assertEqual(S, net.stoich_matrix)
        Y = sp.Matrix([[1, 0, 0, 0], [0, 2, 1, 0], [0, 0, 0, 1], [0, 0, 0, 1]])
        self.assertEqual(Y, net.complex_matrix)
        Ia = sp.Matrix([[-1,  0,  0], [ 1,  0,  0], [ 0, -1,  1], [ 0,  1, -1]])
        self.assertEqual(Ia, net.incidence_matrix)
        k_r0, k_r1, k_r1_rev = sp.symbols("k_r0, k_r1, k_r1_rev")
        self.assertEqual((k_r0, k_r1, k_r1_rev), net.kinetic_params)
        L = sp.Matrix([[k_r0, 0, 0, 0], [-k_r0, 0, 0, 0], [0, 0, k_r1, -k_r1_rev], [0, 0, -k_r1, k_r1_rev]])
        self.assertEqual(L, net.laplacian)
        self.assertEqual((), net.removed_species)

        net = from_react_strings(["A -> 2B", "B <-> C + D"])
        self.assertEqual(None, net.model)
        net.update_model()
        self.assertNotEqual(None, net.model)
        self.assertEqual(['A', 'B', 'C', 'D'], [net.model.getSpecies(s).getId() for s in range(net.model.getNumSpecies())])
        self.assertEqual(("A", "B", "C", "D"), net.species)
        self.assertEqual(4, net.n_species)
        self.assertEqual(4, net.n_complexes)
        self.assertEqual(3, net.n_reactions)
        self.assertEqual(("r0", "r1", "r1_rev"), net.reactionids)
        self.assertEqual(sorted(["A", "2B", "B", "C + D"]), sorted(list(map(str, net.complexes))))
        self.assertEqual(("r0: A ->(k_r0) 2B", "r1: B ->(k_r1) C + D", "r1_rev: C + D ->(k_r1_rev) B"), tuple(map(str, net.reactions)))
        self.assertEqual(sp.Matrix([sympify("k_r0*A"), sympify("k_r1*B"), sympify("k_r1_rev*C*D")]), net.rates)
        S = sp.Matrix([[-1,  0,  0], [ 2, -1,  1], [ 0,  1, -1], [ 0,  1, -1]])
        self.assertEqual(S, net.stoich_matrix)
        Y = sp.Matrix([[1, 0, 0, 0], [0, 2, 1, 0], [0, 0, 0, 1], [0, 0, 0, 1]])
        self.assertEqual(Y, net.complex_matrix)
        Ia = sp.Matrix([[-1,  0,  0], [ 1,  0,  0], [ 0, -1,  1], [ 0,  1, -1]])
        self.assertEqual(Ia, net.incidence_matrix)
        k_r0, k_r1, k_r1_rev = sp.symbols("k_r0, k_r1, k_r1_rev")
        L = sp.Matrix([[k_r0, 0, 0, 0], [-k_r0, 0, 0, 0], [0, 0, k_r1, -k_r1_rev], [0, 0, -k_r1, k_r1_rev]])
        self.assertEqual(L, net.laplacian)
        self.assertEqual((), net.removed_species)

        net.reactions = parse_reactions(["r_id_1: E -> F", "r_id_2: F + A -> 3G"])
        self.assertNotEqual(None, net.model)
        self.assertEqual(('A', 'E', 'F', 'G'), net.species)
        # check that model is updated
        self.assertEqual(['A', 'E', 'F', 'G'], [net.model.getSpecies(s).getId() for s in range(net.model.getNumSpecies())])
        self.assertEqual(4, net.n_species)
        self.assertEqual(4, net.n_complexes)
        self.assertEqual(2, net.n_reactions)
        self.assertEqual(("r_id_1", "r_id_2"), net.reactionids)
        self.assertEqual(sorted(["A + F", "E", "F", "3G"]), sorted(list(map(str, net.complexes))))
        self.assertEqual(("r_id_1: E ->(k_r_id_1) F", "r_id_2: A + F ->(k_r_id_2) 3G"), tuple(map(str, net.reactions)))
        k_r_id_1, k_r_id_2 = sp.symbols("k_r_id_1, k_r_id_2")
        self.assertEqual((k_r_id_1, k_r_id_2), net.kinetic_params)
        self.assertEqual(sp.Matrix([sympify("k_r_id_1*E"), sympify("k_r_id_2*A*F")]), net.rates)
        S = sp.Matrix([[ 0, -1], [-1,  0], [ 1, -1], [ 0,  3]])
        self.assertEqual(S, net.stoich_matrix)
        Y = sp.Matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 1, 0], [0, 0, 0, 3]])
        self.assertEqual(Y, net.complex_matrix)
        Ia = sp.Matrix([[-1,  0], [ 1,  0], [ 0, -1], [ 0,  1]])
        self.assertEqual(Ia, net.incidence_matrix)
        L = sp.Matrix([[k_r_id_1, 0, 0, 0], [-k_r_id_1, 0, 0, 0], [0, 0, k_r_id_2, 0], [0, 0, -k_r_id_2, 0]])
        self.assertEqual(L, net.laplacian)
        self.assertEqual((), net.removed_species)


    def test_save_sbml(self):
        crn = from_react_strings(["r_123: A -> B + C", "2A ->(0.0012345) E", "r_in: -> E"])
        print crn.model
        output = path.join(input_sbml, "test_save_sbml.xml")
        compare = path.join(input_sbml, "test_save_sbml_check.xml")
        crn.save_sbml(output)
        self.assertTrue(cmp(output, compare))


    def test_save_react_file(self):
        crn = from_react_strings(["r_123: A -> B + C", "2A ->(0.0012345) E", "r_in: -> E"])
        output1 = path.join(input_reactions, "test_save_reacts_1")
        output2 = path.join(input_reactions, "test_save_reacts_2")
        compare1 = path.join(input_reactions, "test_save_reacts_check_1")
        compare2 = path.join(input_reactions, "test_save_reacts_check_2")
        crn.save_reaction_file(output1)
        self.assertTrue(cmp(output1, compare1))
        crn.save_reaction_file(output2, rate = True, precision = 2)
        self.assertTrue(cmp(output2, compare2))


    def test_update(self):
        crn = from_react_strings(['a ->(k1) b', 'b (k_2)<->(k2) c + d', 'c ->(k3) d'])
        self.assertEqual(('a', 'b', 'c', 'd'), crn.species)
        self.assertEqual(['a', 'b', 'c + d', 'c', 'd'], [str(c) for c in crn.complexes])
        S = sp.Matrix([[-1, 0, 0, 0], [1, -1, 1, 0], [0, 1, -1, -1], [0, 1, -1, 1]])
        k1, k2, k_2, k3, k4, k_4 = sp.symbols('k1 k2 k_2 k3 k4 k_4')
        L = sp.Matrix([[k1, 0, 0, 0, 0], [-k1, k2, -k_2, 0, 0], [0, -k2, k_2, 0, 0],
                       [0, 0, 0, k3, 0], [0, 0, 0, -k3, 0]])
        self.assertEqual(S, crn.stoich_matrix)
        self.assertEqual(L, crn.laplacian)

        crn.reactions = parse_reactions(['a ->(k1) b', 'c + e (k_4)<->(k4) f'])
        self.assertEqual(('a', 'b', 'c', 'e', 'f'), crn.species)
        self.assertEqual(None, crn.model)
        crn.update_model()
        self.assertNotEqual(None, crn.model)
        self.assertEqual(['a', 'b', 'c', 'e', 'f'], [crn.model.getSpecies(s).getId() for s in range(crn.model.getNumSpecies())])
        self.assertEqual(3, crn.model.getNumReactions())
        self.assertEqual(['a', 'b', 'c + e', 'f'], [str(c) for c in crn.complexes])
        S = sp.Matrix([[-1, 0, 0], [1, 0, 0], [0, -1, 1], [0, -1, 1], [0, 1, -1]])
        L = sp.Matrix([[k1, 0, 0, 0], [-k1, 0, 0, 0], [0, 0, k4, -k_4], [0, 0, -k4, k_4]])
        self.assertEqual(S, crn.stoich_matrix)
        self.assertEqual(L, crn.laplacian)


    def test_is_ma(self):
        crn = from_react_strings(['a ->(k1) b', '2a ->(k2) '])
        self.assertTrue(crn.is_ma)
        crn = from_react_strings(['a ->(k1*a) b'])
        self.assertFalse(crn.is_ma)
        crn = from_react_strings(['a ->(k1*a/a) b', '-> d'])
        self.assertTrue(crn.is_ma)
        crn = from_react_strings(['s ->(k1/(s+1)) p'])
        self.assertFalse(crn.is_ma)


    def test_intermediates(self):
        reacts = ['a -> y', 'y + b -> c']
        crn = from_react_strings(reacts)
        self.assertTrue(crn.is_intermediate_species('y'))
        self.assertTrue(crn.has_linear_equation('y'))
        self.assertTrue('y' in crn.intermediate_species)
        self.assertTrue('y' in crn.intermediate_stoich_1_species)

        reacts = ['a -> y', 'y + b ->(k*y) c']
        crn = from_react_strings(reacts)
        self.assertTrue(crn.is_intermediate_species('y'))
        self.assertFalse(crn.has_linear_equation('y'))
        self.assertTrue('y' in crn.intermediate_species)
        self.assertTrue('y' in crn.intermediate_stoich_1_species)

        reacts = ['a -> y', 'y + b -> y + c']
        crn = from_react_strings(reacts)
        self.assertFalse(crn.is_intermediate_species('y'))
        self.assertFalse(crn.has_linear_equation('y'))
        self.assertFalse('y' in crn.intermediate_species)
        self.assertFalse('y' in crn.intermediate_stoich_1_species)

        reacts = ['a -> 2y', '2y + b -> c']
        crn = from_react_strings(reacts)
        self.assertTrue(crn.is_intermediate_species('y'))
        self.assertFalse(crn.has_linear_equation('y'))
        self.assertTrue('y' in crn.intermediate_species)
        self.assertFalse('y' in crn.intermediate_stoich_1_species)


    def test_source_sink(self):
        reacts = ['a -> y', 'b <-> y', 'd -> c', 'y + d -> c', 'd -> e + d']
        crn = from_react_strings(reacts)
        self.assertEqual(['a', 'd'], crn.source_species)
        self.assertEqual(['c', 'e'], crn.sink_species)
        self.assertEqual(['a', 'd', 'd + y'], list(map(str, crn.source_complexes)))
        self.assertEqual(['c', 'd + e'], list(map(str, crn.sink_complexes)))
        self.assertEqual(['y', 'b'], list(map(str, crn.intermediate_complexes)))
        self.assertEqual(['b'], list(map(str, crn.intermediate_simple_complexes)))


    def test_conn_comp(self):
        reacts = ['a -> b', 'b -> c', 'c -> a', 'd <-> e', 'e -> f']
        crn = from_react_strings(reacts)
        self.assertEqual(2, crn.weak_conn_components()[0])
        self.assertEqual(3, crn.strong_conn_components()[0])
        self.assertEqual(2, len(crn.strong_terminal_conn_components()[0]))


    def test_detach(self):
        crn = from_react_file(path.join(input_reactions, "test_detach"))
        net2 = crn.detach_network([crn.complexes.index(to_complex('D'))])
        self.assertEqual(set(['P', 'F', 'Q', 'D']), set(net2.species))
        net2.qss(cons_law = ('F', ConsLaw('F + D', 'Ftot')))
        for s, f in net2.removed_species: print(s + " = " + str(f))
        net2.inspect(True, print_matrices = True)

        crn = from_react_file(path.join(input_reactions, "test_detach"))
        crn.qss(cons_law = ('F', ConsLaw('F + D', 'Ftot')))
        crn.inspect(True, print_matrices = True)


    def test_dynEq(self):
        crn1 = from_react_file(path.join(input_reactions, "dynEq/net1"))
        crn2 = from_react_file(path.join(input_reactions, "dynEq/net2"))
        crn1.print_equations()
        crn2.print_equations()
        assert crn1.is_dyn_eq(crn2)
        assert crn2.is_dyn_eq(crn1)


    def test_emul(self):
        crn1 = from_react_file(path.join(input_reactions, "emul/emul_not_isom_1"))
        crn2 = from_react_file(path.join(input_reactions, "emul/emul_not_isom_2"))

        crn1.print_equations()
        crn2.print_equations()

        self.assertRaises(ValueError, crn2.is_emul, crn1)
        self.assertRaises(ValueError, crn2.is_emul, crn1, morphism = {"A": "A"})
        self.assertFalse(crn2.is_emul(crn1, morphism = {"A": "A", "B": "B", "C": "A"}))
        self.assertTrue(crn2.is_emul(crn1, morphism = {"A": "A", "B": "B", "C": "B"}))
        self.assertFalse(crn1.is_dyn_eq(crn2))


    def test_dynEq2(self):
        reacts1 = ["A ->(k1) B", "B ->(k_1) A", "B ->(k2) C", "C ->(k_2) B"]
        crn1 = from_reacts(parse_reactions(reacts1))
        reacts2 = ["->(k_1*(A + B + C)) A", "A ->(k1 + k_1 + k_1 * C/A)", "->(k2*(A + B + C)) C", "C ->(k2 + k_2 + k2 * A/C)"]
        crn2 = from_reacts(parse_reactions(reacts2))
        crn1.print_equations()
        crn2.print_equations()
        self.assertTrue(crn1.is_emul(crn2))
        self.assertFalse(crn1.is_dyn_eq(crn2))


    def test_derivative(self):
        reacts = ['a -> b + c', '2c -> a + d']
        crn = from_react_strings(reacts)
        crn.print_equations()
        print
        self.assertRaises(ValueError, crn.derivative, 'b + 3')
        self.assertRaises(ValueError, crn.derivative, 'a**2 + b')
        self.assertRaises(ValueError, crn.derivative, '2 a + f')
        self.assertEqual(0, (sp.sympify('k_r0*a + k_r1*c**2') - crn.derivative('a + 2b')).simplify())


    def test_im(self):
        crn = from_react_file(path.join(input_reactions, "dsr-graph/pos_loops_main"))
        g_1_1, g_2_2, g_2_3, g_3_3, g_3_4 = sp.symbols('g_1_1 g_2_2 g_2_3 g_3_3 g_3_4')
        I = sp.Matrix([[g_1_1, 0, 0, 0], [0, g_2_2, g_2_3, 0], [0, 0, g_3_3, g_3_4]])
        self.assertEqual(I, crn.influence_matrix(state = {'x1': 2}))
        crn.dsr_graph_adj()
        crn = from_react_file(path.join(input_reactions, "dsr-graph/pos_loops_main_v2"))
        I = sp.Matrix([[g_1_1, 0, 0, 0], [0, g_2_2, g_2_3, 0], [0, 0, 0, g_3_4]])
        self.assertEqual(I, crn.influence_matrix(check = True))
        crn.dsr_graph_adj()
        a1_1, a2_2, a2_3, a3_3, a3_4 = sp.symbols('a1_1 a2_2 a2_3 a3_3 a3_4')
        I = sp.Matrix([[a1_1, 0, 0, 0], [0, a2_2, a2_3, 0], [0, 0, 0, a3_4]])
        self.assertEqual(I, crn.influence_matrix(var = "a", check = True, params = {'K': 7, 'k1': 3}))
        crn = from_react_strings(["x ->(k/(1-x))"])
        self.assertRaises(ValueError, crn.influence_matrix, check = True)


    def test_fix_ma(self):
        reacts = ["b ->(k1*b) c", "a ->(k2*b**2) c"]
        reacts_fixed = ["r0: 2b ->(k1) b + c", "r1: a + 2b ->(k2) 2b + c"]
        crn = from_reacts(parse_reactions(reacts))
        for r in crn.reactions: print(r)
        print
        crn._fix_ma()
        for r in crn.reactions: print(r)
        print
        self.assertEqual([str(r) for r in crn.reactions], reacts_fixed)


    def test_fix_denom(self):
        reacts = ["a + b ->(k1/b) b", "a + 3 b ->(k2/b**2) c + 5b"]
        reacts_fixed = ["r0: a ->(k1) ", "r1: a + b ->(k2) 3b + c"]
        crn = from_reacts(parse_reactions(reacts))
        for r in crn.reactions: print(r)
        print
        crn._fix_denom()
        for r in crn.reactions: print(r)
        print
        self.assertEqual([str(r) for r in crn.reactions], reacts_fixed)


    def test_acr(self):
        crn = from_react_file(path.join(input_reactions, "acr/acr_1"))
        self.assertEqual(['yp'], crn.acr_species())

        crn = from_react_file(path.join(input_reactions, "acr/acr_toy"))
        self.assertEqual(['a'], crn.acr_species())

        crn = from_react_file(path.join(input_reactions, "acr/acr_complex"))
        self.assertEqual([sp.sympify('A*B')], crn.acr_complexes())
        self.assertEqual(['C'], crn.acr_species(subnets = True))

        crn = from_react_file(path.join(input_reactions, "acr/neigenfind_ex1"))
        self.assertEqual(['C'], crn.acr_species(subnets = True))

        crn = from_react_file(path.join(input_reactions, "acr/neigenfind_ex2"))
        self.assertEqual(['C', 'U'], crn.acr_species(subnets = True))


    def test_tree_constants(self):
        crn = from_react_strings(['A1 (k2)<->(k1) A2', 'A2 (k4+k1*k5/k2)<->(k3) A3'])
        tconstants = list(map(sp.sympify, ['k2*k4+k1*k5', 'k1/k2*(k2*k4+k1*k5)', 'k1*k3']))
        tree_consts = crn.tree_constants()
        for i in range(crn.n_complexes):
            self.assertEqual((crn.tree_constant(i) - tconstants[i]).simplify(), 0)
        for i in range(crn.n_complexes):
            self.assertEqual((tree_consts[i] - tconstants[i]).simplify(), 0)


    def test_ems(self):
        reacts = ["a -> b", "b + c -> d", "d <-> e", "e <-> a + c", "e -> b + c"]
        crn = from_react_strings(reacts)
        self.assertRaises(ValueError, crn.is_ss_flux, sp.Matrix([0, 0, 2, 2, 3, 3]))
        self.assertTrue(crn.is_ss_flux(sp.Matrix([0, 0, 2, 2, 3, 3, 0])))
        self.assertFalse(crn.is_ss_flux(sp.Matrix([1, 1, 1, 1, 1, 1, 1])))
        self.assertRaises(ValueError, crn.is_cyclic_ss_flux, sp.Matrix([0, 0, 2, 2, 3, 3]))
        self.assertTrue(crn.is_cyclic_ss_flux(sp.Matrix([0, 0, 2, 2, 3, 3, 0])))
        self.assertFalse(crn.is_cyclic_ss_flux(sp.Matrix([1, 1, 1, 0, 3, 2, 0])))
        self.assertEqual([False, True, True, True], [crn.is_cyclic_ss_flux(v) for v in crn.t_invariants.tolist()])


    def test_set_params(self):
        reacts = ["r1: a ->(k1) b", "r2: b + c ->(k2) d", "d (k_3)<->(k3) a + c"]
        crn = from_react_strings(reacts)
        crn.set_params({'k1': 0.001, 'k2': 1, 'k_3': 'k3**2/d'})
        self.assertEqual(str(crn.reactions[0]), "r1: a ->(1.000e-3) b")
        self.assertEqual(str(crn.reactions[1]), "r2: b + c ->(1) d")
        self.assertEqual(str(crn.reactions[2]), "r0: d ->(k3) a + c")
        self.assertEqual(str(crn.reactions[3]), "r0_rev: a + c ->(k3**2/d) d")
        crn.set_params({'k3': 5e7})
        for r in crn.reactions: print(r)
        self.assertEqual(str(crn.reactions[0]), "r1: a ->(1.000e-3) b")
        self.assertEqual(str(crn.reactions[1]), "r2: b + c ->(1) d")
        self.assertEqual(str(crn.reactions[2]), "r0: d ->(5.000e+7) a + c")
        self.assertEqual(str(crn.reactions[3]), "r0_rev: a + c ->(2.5e+15/d) d")


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCrn)
    unittest.TextTestRunner(verbosity=2).run(suite)
