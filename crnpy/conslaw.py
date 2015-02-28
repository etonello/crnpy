#!/usr/bin/env python

"""Conservation Law class."""

import logging

from .crncomplex import sympify

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


class ConsLaw:
    """Conservation law of the form expression = constant,
    for example e + es + esi = etot.

    The parameters expression and constant are strings or sympy expressions.

    Attributes:
        * expression: the conservation expression converted to sympy expression.
        * constant: the conservation constant converted to sympy expression.

    :Example:

    >>> cons_law = ConsLaw("E + ES + ESI", "Etot")
    >>> cons_law
    E + ES + ESI = Etot
    >>> cons_law.expression.subs("ESI", "ES*I")
    E + ES + I*ES
    >>> cons_law.constant
    Etot

    """
    def __init__(self, expression, constant):
        self.logger = logging.getLogger("crn.conservation_law")
        self.logger.info("Creating an instance of conservation law.")
        self.expression = sympify(expression)
        self.constant = sympify(constant)
        self.species = self.expression.as_coefficients_dict()

    def __str__(self):
        return str(self.expression) + " = " + str(self.constant)

    def  __repr__(self):
        return self.__str__()
