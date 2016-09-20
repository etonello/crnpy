#!/usr/bin/env python

"""Tests for Reaction class."""

import unittest

from crnpy.crn import CRN, from_react_file
from crnpy.crncomplex import *
from crnpy.parsereaction import parse_complex

from sympy import SympifyError

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


class TestComplex(unittest.TestCase):
    def test_to_complex(self):
        self.assertEqual(to_complex('2*a + b'), Complex({'a': 2, 'b': 1}))
        self.assertRaises(SympifyError, to_complex, '2a + b')
        self.assertRaises(ValueError, to_complex, '2*a + 3')


    def test_ma(self):
        self.assertEqual(to_complex('2*a + b').ma(), sympify('a**2*b'))
        self.assertEqual(to_complex('a + 3*d + b').ma(), sympify('d**3*a*b'))
        self.assertEqual(parse_complex('2a + b').ma(), sympify('a**2*b'))
        self.assertEqual(parse_complex('a + 3d + b').ma(), sympify('d**3*a*b'))
        self.assertEqual(parse_complex('').ma(), sympify('1'))


    def test_to_vector(self):
        self.assertEqual(to_complex('2*a + b').to_vector(['s', 'a', 'b', 'c']),
                         sp.Matrix([0, 2, 1, 0]))


    def test_symp(self):
        self.assertEqual(Complex({'a': 1, 'b': 2, 's': 1}).symp(), sp.sympify('a + 2*b + s'))


    def test_times(self):
        self.assertEqual(Complex({'a': 1, 'b': 2, 's': 1}).times(3),
                         Complex({'a': 3, 'b': 6, 's': 3}))
        self.assertEqual(Complex({'a': 1, 'b': 2, 's': 1}).times(0),
                         Complex({}))
        self.assertEqual(Complex({'a': 1, 'b': 2, 's': 1}).times(-1),
                         Complex({'a': -1, 'b': -2, 's': -1}))
        self.assertRaises(ValueError, Complex({'a': 1, 'b': 2}).times, '2*a')
        self.assertRaises(ValueError, Complex({'a': 1, 'b': 2}).times, '2')
        self.assertRaises(ValueError, Complex({'a': 1, 'b': 2}).times, sp.sympify('2*a'))
        self.assertRaises(ValueError, Complex({'a': 1, 'b': 2}).times, sp.sympify('2'))


    def test_order(self):
        self.assertTrue(Complex({'a': 1, 'b': 1}) <= Complex({'a': 1, 'b': 2, 's': 1}))
        self.assertTrue(Complex({'a': 1, 'b': 1}) < Complex({'a': 1, 'b': 2, 's': 1}))
        self.assertFalse(Complex({'a': 1, 'b': 1}) >= Complex({'a': 1, 'b': 2, 's': 1}))
        self.assertFalse(Complex({'a': 1, 'b': 1}) > Complex({'a': 1, 'b': 2, 's': 1}))
        self.assertFalse(Complex({'a': 1, 'b': 1}) > Complex({'a': 1, 'c': 1}))
        self.assertTrue(Complex({'a': 1, 'b': 1}) > Complex({'a': 1}))
        self.assertTrue(Complex({'a': 1, 'b': 1}) < Complex({'a': 1, 'b': 1, 'c': 3}))
        self.assertTrue(Complex({'a': 1, 'b': 1}) <= Complex({'a': 1, 'b': 1, 'c': 3}))
        self.assertFalse(Complex({'a': 1, 'b': 1}) <= Complex({'a': 1, 'c': 3}))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestComplex)
    unittest.TextTestRunner(verbosity=2).run(suite)
