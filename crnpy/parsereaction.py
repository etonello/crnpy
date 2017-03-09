#!/usr/bin/env python

"""Parser for reaction strings and files."""

import libsbml
import re
import sympy as sp
from sympy.abc import _clash

from .crncomplex import Complex
from .reaction import Reaction

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


def parse_reaction(r):
    """Parse a (potentially reversible) reaction written as

    reactionid: i1 R1 + ... + in Rn (k\_)<->(k) j1 P1 + ... + jm Pm

    where (k\_)< is optional.
    Return a two-element tuple with two reactions, if "<" is included
    (the reaction is reversible),
    or one reaction and None, if it is not reversible.
    k, k\_ can be any string representing the rates. They are optional. If present,
    they must be enclosed in parenthesis.
    No spaces are allowed between the parenthesis and <->.
    reactionid is optional. If the reaction is reversible, "_rev"
    is added to the reactionid for the reverse reaction.
    Everything after a # sign is ignored.

    :type r: string
    :param r: string of the form "reactionid: i1 R1 + ... + in Rn (k\_)<->(k) j1 P1 + ... + jm Pm".
    :rtype: (Reaction, Reaction/None).
    """
    # Everything after a # sign is ignored.
    reaction = r.split("#")[0]

    # Check for reaction id.
    colon = reaction.split(":")
    if len(colon) > 2:
        raise ValueError("Unrecognised reaction. More then one colon in reaction definition.")
    if len(colon) == 2:
        reactionid = colon[0].strip()
        reactionBody = colon[1]
    else:
        reactionid = None
        reactionBody = colon[0]

    if reactionBody != "":
        pattern = re.compile("(.*?)(?:\((.*)\))?(<)?\->(?:\((.*)\))?(.*)")
        m = pattern.match(reactionBody)
        if m is None:
            raise ValueError("Unrecognised reaction.")
        else:
            reacts, k_, inv, k, prods = m.groups()

        reactants = parse_complex(reacts)
        products = parse_complex(prods)

        if inv == "<":
            if reactionid != None: reactionidRev = reactionid + "_rev"
            else: reactionidRev = None
            return (Reaction(reactionid, reactants, products, parse_expr(k)), \
                    Reaction(reactionidRev, Complex(products), Complex(reactants), parse_expr(k_)))
        else:
            return (Reaction(reactionid, reactants, products, parse_expr(k)), None)


def _valid_species(cexpr):
    try:
        sp.Symbol(cexpr)
    except:
        raise ValueError("Could not parse complex {}".format(cexpr))


def parse_complex(complex_string):
    """Parse a string representing a complex.
    The string must be of the form "n1 s1 + n2 s2 + ...",
    where the ni are the integer stoichiometric coefficients.
    Stoichiometric coefficients equal to 1 can be omitted.

    :type complex_string: string
    :param complex_string: string of the form "n1 s1 + n2 s2 + ...".
    :rtype: Complex.
    """
    complex_string = complex_string.replace(" ", "").split("+")

    pattern = re.compile("(\d*)(?:\*)?(.*)")

    parsedComplex = {}

    for c in complex_string:
        m = pattern.match(c)
        if m is None: raise ValueError("Unrecognised complex.")
        m = m.groups()
        if m[1] != '': _valid_species(m[1])
        if m[0] == '' and m[1] != '':
            parsedComplex[m[1]] = 1
        else:
            if m[0] != '':
                parsedComplex[m[1]] = int(m[0])

    return Complex(parsedComplex)


def param_to_rate(reaction):
    """Multiplies the rate by the reactant mass-action monomial.

    :type reaction: Reaction
    :rtype: Reaction.
    """
    return Reaction(reaction.reactionid, \
                    reaction.reactant, \
                    reaction.product, \
                    reaction.rate * reaction.reactant.ma())


def _read_reactions(reacts):
    """Parse each reaction into forward
    and backward reactions."""
    rs = []
    for row in reacts:
        row = row.strip()
        if row:
            r = parse_reaction(row)
            if r:
                rs.append(r)
    return rs


def parse_reaction_file(filename, rate = False):
    """Read a reaction file, and populate
    reaction id and kineticParam if empty.
    If rate is False, the expression enclosed in parenthesis
    is interpreted as rate / reactant mass-action monomial,
    otherwise as the full rate.

    :type filename: string
    :type rate: boolean
    :param filename: path to file.
    :param rate: True if the expressions in parenthesis are rates.
    :rtype: list of Reactions.
    """
    with open(filename) as f:
        reactions = parse_reactions(f.readlines(), rate)
    return reactions


def parse_reactions(rs, rate = False):
    """Parse a list of reaction strings.
    If rate is False, the expression enclosed in parenthesis
    is interpreted as rate / reactant mass-action monomial,
    otherwise as the full rate.

    :type rs: list of strings
    :param rs: strings representing the reactions.
    :rtype: list of Reactions.
    """
    if not isinstance(rs, list):
        raise ValueError("Required list of strings")
    if rate:
        return list(map(add_kinetic_param, add_reaction_id(_read_reactions(rs))))
    else:
        reacts = _read_reactions(rs)
        return list(map(param_to_rate, list(map(add_kinetic_param, add_reaction_id(reacts)))))


def add_reaction_id(reactions):
    """ Add a reactionid, if missing, of the form
    'ri', where i is the index of the reaction, in the sublist of reactions without reactionid.
    If there is a reverse reaction, its id is set to ri_rev.

    :type reactions: list of Reactions.
    :rtype: list of Reactions.
    """

    rids = [r[0].reactionid for r in reactions if r[0].reactionid]
    if len(rids) > 0 and len(rids) != len(list(set(rids))):
        raise ValueError("Non-unique reaction ids.")

    get_id = ("r" + str(j) for j in range(len(reactions))
              if "r" + str(j) not in rids)

    newreactions = []
    for i in range(len(reactions)):
        reaction, rev_reaction = reactions[i]
        if not reaction.reactionid:
            reaction._reactionid = next(get_id)
            newreactions.append(reaction)
            if rev_reaction:
                rev_reaction._reactionid = reaction.reactionid + "_rev"
                newreactions.append(rev_reaction)
        else:
            newreactions.append(reaction)
            if rev_reaction: newreactions.append(rev_reaction)
    return newreactions


def add_kinetic_param(reaction):
    """Add a kinetic param of the form 'k_reactionid', if missing.

    :type reaction: Reaction.
    :rtype: Reaction.
    """
    if not reaction._rate: reaction._rate = sp.Symbol("k_" + reaction.reactionid)
    return reaction


def ast_to_sympy_expr(math, debug = False):
    """Create a sympy expression from a libsbml AST.

    :type math: libsbml.ASTNode
    :rtype: sympy expression
    """
    t = math.getType()
    if debug:
        print("Node of type {}, {}".format(dtypes[t], math.getName()))
    if t == libsbml.AST_NAME:
        if debug:
            print("Node name: {}".format(math.getName()))
        return sp.Symbol(math.getName())

    if t == libsbml.AST_NAME_TIME:
        return sp.Symbol('time')

    if t == libsbml.AST_CONSTANT_PI:
        # using same symbol as in SBML formula, rather than sp.pi
        return sp.Symbol("Pi")
    if t == libsbml.AST_CONSTANT_E:
        return sp.E

    if math.isInteger():
        return sp.Integer(math.getInteger())
    if math.isRational():
        return sp.Rational(math.getNumerator(), math.getDenominator())
    if math.isReal():
        return sp.Float(math.getReal())

    nc = math.getNumChildren()
    children = [ast_to_sympy_expr(math.getChild(c)) for c in range(nc)]

    if t == libsbml.AST_PLUS:
        return sp.Add(*children)

    if t == libsbml.AST_TIMES:
        return sp.Mul(*children)

    if t == libsbml.AST_POWER or t == libsbml.AST_FUNCTION_POWER:
        if len(children) != 2:
            raise ValueError("Error parsing SBML kineticLaw: {}", t)
        return children[0]**children[1]

    if t == libsbml.AST_MINUS:
        if len(children) != 1 and len(children) != 2:
            raise ValueError("Error parsing SBML kineticLaw: {}", t)
        if len(children) == 1:
            return - children[0]
        else:
            return children[0] - children[1]

    if t == libsbml.AST_DIVIDE:
        if len(children) != 2:
            raise ValueError("Error parsing SBML kineticLaw: {}", t)
        return children[0] / children[1]

    if t == libsbml.AST_FUNCTION_EXP:
        if len(children) != 1:
            raise ValueError("Error parsing SBML kineticLaw: {}", t)
        return sp.E**children[0]

    if t == libsbml.AST_FUNCTION_DELAY:
        return sp.Function('delay')(*children)

    if math.isFunction():
        return sp.Function(math.getName())(*children)

    raise NotImplementedError("Type {} of AST Node not supported.", t)


def parse_expr(s):
    """Try to convert a string to a sympy expression,
    reading the string using the SBML formula parser first,
    and converting the libSBML AST to a sympy expression.
    """
    if not s:
        return None
    s = s.replace('**', '^')
    ast_tree = libsbml.parseL3Formula(s)
    if ast_tree:
        return ast_to_sympy_expr(ast_tree)
    else:
        raise ValueError("Could not parse expression.")


def flux_value(math):
    if math.getType() == libsbml.AST_NAME and \
       math.getName()[-10:] == "FLUX_VALUE":
        return math.getName()
