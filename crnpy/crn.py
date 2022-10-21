#!/usr/bin/env python

"""CRN class."""

from collections import defaultdict
import itertools
from libsbml import writeSBMLToFile, formulaToL3String, SBMLDocument
import logging
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import sympy as sp
import warnings
import sympy as sp
from scipy import integrate

from .conslaw import ConsLaw
from .createmodel import model_from_reacts, model_from_sbml, replace_reacts
from .matrixfunctions import negative, sdiag, print_matrix, _pos_dependent, _pos_generators
from .crncomplex import Complex
from .parsereaction import parse_reaction_file, parse_complex, parse_reactions, ast_to_sympy_expr, flux_value
from .reaction import Reaction, _split_reaction_monom, merge_reactions, _same_denom

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


class CRN(object):
    """A chemical reaction network.

    Keyword argument:
        * model: SBML model.

    :Example:

    >>> net = CRN(model_from_sbml("enzyme.xml"))

    or equivalently

    >>> net = from_sbml("enzyme.xml")

    In alternative, a CRN can be created from a list of reactions or strings representing reactions:

    >>> reacts = [Reaction("r1", Complex(a = 1), Complex(b = 1), "k1*a"), Reaction("r2", Complex(b = 2), Complex(), "k2*b**2")]
    >>> net = from_reacts(reacts)

    or equivalently

    >>> net1 = from_react_strings(["r1: a ->(k1) b", "r2: 2b ->(k2) "])

    One can also create an empty CRN and update the reactions:

    >>> from crnpy.parsereaction import parse_reactions
    >>> net2 = CRN()
    >>> net2.species
    ()
    >>> net2.complexes
    ()
    >>> net2.reactions = parse_reactions(["r1: a ->(k1) b", "r2: 2b ->(k2) "])
    >>> net2.species
    ('a', 'b')
    >>> net2.complexes
    (a, b, 2b, )
    """

    def __init__(self, model = None):
        self.logger = logging.getLogger("crnpy.crn")
        self.logger.info("Creating an instance of crn.")

        if model:
            self._model, self._document = model

            # species
            self._species_from_sbml()
            self._removed_species = []

            # reactions
            self._get_reactions()
            self._populate()

        else:
            self._model, self._document = None, None
            self._removed_species = []
            self.reactions = []


    @classmethod
    def from_reacts(cls, reacts):
        """Create a reaction network from a list of reactions."""
        crn = cls()
        crn.reactions = reacts
        return crn


    @property
    def model(self):
        """SBML model of the chemical reaction network.

        If the reaction network is created using a list of reactions,
        the SBML model is not created. It is created if save_sbml is called,
        or if update_model() is called.
        """
        return self._model


    @property
    def document(self):
        """SBML document.

        If the reaction network is created using a list of reactions,
        the SBML document is not created. It is created if save_sbml is called,
        or if update_model() is called.
        """
        return self._document


    @property
    def species(self):
        """Tuple of the network species.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 ->(k2) 2 A3"])
        >>> net.species
        ('A1', 'A2', 'A3')

        The tuple of species is read-only, and can change if the reactions are updated,
        or if a reduction method is applied to eliminate some species.

        :rtype: tuple of strings.
        """
        return tuple(self._species)


    @property
    def removed_species(self):
        """Pairs (species, expression) for species that have been eliminated.

        :rtype: list of pairs (string, sympy expression).
        """
        return tuple(self._removed_species)


    @property
    def complexes(self):
        """Tuple of the network complexes.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 ->(k2) 2 A3"])
        >>> net.complexes
        (A1, A2 + A3, A2, 2A3)

        The tuple of complexes is read-only, and can change if the reactions are updated,
        or if a reduction method is applied to eliminate some species.

        :rtype: tuple of Complex objects.
        """
        return tuple(self._complexes)


    @property
    def reactions(self):
        """Tuple of the network reactions.

        :setter: Sets the reactions, updating the CRN model and document if they exist.
        :rtype: tuple of Reaction objects.
        """
        return tuple(self._reactions)
    @reactions.setter
    def reactions(self, value):
        self._species = sorted(list(set([s for r in value for s in r.reactant] +
                                       [s for r in value for s in r.product])))
        self._n_species = len(self._species)
        self._reactions = value
        self._populate()
        self.update_model(if_exists = True)


    @property
    def n_species(self):
        return self._n_species


    @property
    def n_complexes(self):
        return self._n_complexes


    @property
    def n_reactions(self):
        return self._n_reactions


    @property
    def reactionids(self):
        return tuple(self._reactionids)


    @property
    def rates(self):
        """Rates of the reactions.

        :rtype: matrix of sympy expressions.
        """
        return self._rates


    @property
    def kinetic_params(self):
        """Kinetic parameters of the reactions.

        :rtype: tuple of sympy expressions.
        """
        return tuple(self._kinetic_params)


    def _species_from_sbml(self):
        """Extract species from SBML model."""
        # species
        self._species = [self.model.getSpecies(s).getId() for s in range(self.model.getNumSpecies())]
        self._species = sorted(self.species)
        self._n_species = len(self.species)


    def _get_reactions(self):
        """Extract reactions from SBML model."""
        nr = self.model.getNumReactions()
        reactions = []

        for r in range(nr):
            reaction = self.model.getReaction(r)
            reactionid = reaction.getId()

            reactant = Complex(dict((c.getSpecies(), 1 if np.isnan(c.getStoichiometry()) else int(c.getStoichiometry())) \
                               for c in reaction.getListOfReactants()))
            product = Complex(dict((c.getSpecies(), 1 if np.isnan(c.getStoichiometry()) else int(c.getStoichiometry())) \
                              for c in reaction.getListOfProducts()))

            # remove species with stoichiometric coefficient equal to 0
            reactant = Complex(dict((s, sc) for s, sc in reactant.items() if sc != 0))
            product = Complex(dict((s, sc) for s, sc in product.items() if sc != 0))

            math = reaction.getKineticLaw().getMath()
            if math:
                # Special case for FLUX_VALUE
                fv = flux_value(math)
                if flux_value(math):
                    rate = sp.Symbol(fv) * reactant.ma()
                    if reaction.getReversible():
                        raterev = sp.Symbol(fv + "_rev") * product.ma()
                else:
                    kineticlaw = ast_to_sympy_expr(math)

                    if reaction.getReversible():
                        # if reaction is reversible, we need to split the kineticLaw
                        # into forward and backward formulas.
                        # We expand the kineticLaw and assume that the components relating
                        # to the inverse reaction are those with a minus sign in front.
                        numer, denom = kineticlaw.as_numer_denom()
                        negative = numer.expand().coeff(-1)
                        rate = ((numer + negative) / denom).factor()
                        raterev = (negative / denom).factor()
                    else:
                        rate = kineticlaw
            else:
                param = "k_" + reaction.reactionid
                rate = reactant.ma() * sp.Symbol(param)
                if reaction.getReversible():
                    raterev = product.ma() * sp.Symbol(param + "_rev")

            reactions.append(Reaction(reactionid, reactant, product, rate))
            if reaction.getReversible():
                revid = reactionid + "_rev"
                if not raterev:
                    raterev = product.ma() * sp.Symbol("k_" + revid)
                reactions.append(Reaction(revid, product, reactant, raterev))
        self._reactions = reactions


    def _populate(self):
        """Create crn attributes from reactions."""
        self._n_reactions = len(self.reactions)
        self._reactionids = [r.reactionid for r in self.reactions]
        self._rates = sp.Matrix([r.rate for r in self.reactions])
        self._kinetic_params = [r.kinetic_param for r in self.reactions]

        self._n_species = len(self.species)

        # complexes and matrices
        self._complexes = []
        nc = -1
        incidence = {}
        cmatrix = {}
        for nr in range(self.n_reactions):
            r = self.reactions[nr]

            if r.reactant in self._complexes:
                indc = self._complexes.index(r.reactant)
            else:
                nc = nc + 1
                indc = nc
                self._complexes.append(r.reactant)
                for s in r.reactant: cmatrix[(self.species.index(s), nc)] = r.reactant[s]
            incidence[(indc, nr)] = -1

            if r.product in self._complexes:
                indc = self._complexes.index(r.product)
            else:
                nc = nc + 1
                indc = nc
                self._complexes.append(r.product)
                for s in r.product: cmatrix[(self.species.index(s), nc)] = r.product[s]
            incidence[(indc, nr)] = 1

            if r.reactant == r.product:
                incidence[(indc, nr)] = 0
        self._n_complexes = len(self._complexes)

        # Create the matrix of complex coefficients. s x c
        self._complex_matrix = sp.SparseMatrix(self.n_species, self.n_complexes, cmatrix)

        # Create the incidence matrix. c x r
        self._incidence_matrix = sp.SparseMatrix(self.n_complexes, self.n_reactions, incidence)


    def update_model(self, if_exists = False):
        """Update the SBML model and document or create them if they do not exist.
        If if_exists is set to True, update only if the model already exists."""
        if (not if_exists) and (not self.model):
            document = SBMLDocument(3, 1)
            self._model, self._document = document.createModel(), document
        if self.model:
            self._model, self._document, self._species = \
                replace_reacts(self.model, self.document, self.reactions)


    def set_params(self, params_dict):
        """Replace the parameters used in the reaction rates
        with the values in dictionary *params_dict*.

        *params_dict* is a dictionary with keys the parameters to replace, and
        values the sympy expressions or numeric values that replace them.

        In the following example we set all the kinetic parameters to 0.001:

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 ->(k2) 2 A3"])
        >>> net.set_params(dict((k, 0.001) for k in net.kinetic_params))
        >>> net.reactions
        (r0: A1 ->(1.000e-3) A2 + A3, r1: A2 ->(1.000e-3) 2A3)

        """
        self.reactions = [Reaction(r.reactionid,
                                   r.reactant,
                                   r.product,
                                   r.rate.subs(params_dict)) for r in self.reactions]
        self.update_model(if_exists = True)


    ### Matrices ####

    @property
    def complex_matrix(self):
        """Complex matrix (usually denoted with Y), i.e. the matrix with dimension
        number of species times number of complexes, with element Yij given by
        the stoichiometric coefficient of the i-th species in the j-th complex.

        :rtype: sympy SparseMatrix.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 ->(k2) 2 A3"])
        >>> net.complex_matrix
        Matrix([
        [1, 0, 0, 0],
        [0, 1, 1, 0],
        [0, 1, 0, 2]])

        """
        # ns x nc matrix
        return self._complex_matrix


    @property
    def incidence_matrix(self):
        """Incidence matrix.

        Sometimes denoted as Ia, it is the matrix with dimension
        number of complexes times number of reactions,
        with element at position ij given by
        -1 if the i-th complex is the reactant of the j-th reaction,
        +1 if the i-th complex is the product of the j-th reaction,
        and 0 otherwise.

        :rtype: sympy SparseMatrix.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 ->(k2) 2 A3"])
        >>> net.incidence_matrix
        Matrix([
        [-1,  0],
        [ 1,  0],
        [ 0, -1],
        [ 0,  1]])

        """

        # nc x nr matrix
        return self._incidence_matrix


    @property
    def stoich_matrix(self):
        """Stoichiometric matrix of the reaction network.

        Matrix with dimension number of species times number of reactions,
        with element at position ij given by
        the difference between the stoichiometric coefficient of the i-th species
        in the product of the j-th reaction and
        the stoichiometric coefficient of the i-th species
        in the reactant of the j-th reaction.

        :rtype: sympy SparseMatrix.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 ->(k2) 2 A3"])
        >>> net.stoich_matrix
        Matrix([
        [-1,  0],
        [ 1, -1],
        [ 1,  2]])

        """
        # s x r matrix
        return sp.SparseMatrix(self.complex_matrix.multiply(self.incidence_matrix))


    @property
    def kinetic_matrix(self):
        """Kinetic matrix.

        Sometimes denoted as Ak, it is the matrix with dimension
        number of complexes times number of complexes,
        with element at position ij given by
        the sum of (generalised) kinetic parameters of reactions
        with reactant the j-th complex and product the i-th complex,
        when i and j are different. The diagonal elements are defined
        so that the columns sum to zero.

        :rtype: sympy SparseMatrix.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 ->(k2) 2 A3"])
        >>> net.kinetic_matrix
        Matrix([
        [-k1, 0,   0, 0],
        [ k1, 0,   0, 0],
        [  0, 0, -k2, 0],
        [  0, 0,  k2, 0]])

        """
        # -Laplacian
        return -self.incidence_matrix.multiply(sdiag(self.kinetic_params).multiply(negative(self.incidence_matrix).T))


    @property
    def laplacian(self):
        """Generalised Laplacian of the graph of complexes.
        It is the negative of the kinetic matrix.

        :rtype: sympy SparseMatrix.
        """
        # -kinetic matrix
        return self.incidence_matrix.multiply(sdiag(self.kinetic_params).multiply(negative(self.incidence_matrix).T))


    def influence_matrix(self, var = None, state = None, check = False, interval = None, params = None):
        """Return the influence matrix of the reaction network.
        This is a matrix of dimension number of species times number of reactions.
        The element at position ij is a variable with a plus in front
        if the rate v_j(x) of reaction j increases with x_i,
        a minus if it decreases, and 0 if it is constant in x_i.
        Every reaction rate v_j is assumed to be strictly monotone in each variable x_i.
        Optional check for monotonicity is available with check = True only for
        functions of one variable.

        :param var: variable name to use. Default is g\_.
        :type var: string
        :param state: dictionary of the point where the derivatives are evaluated. Default is (1, 1, ..., 1).
        :type state: dictionary, key = string, value = numeric
        :param check: check for monotonicity. Only available for unary functions.
        :type check: boolean
        :param interval: interval where monotonicity is checked - default is [0, oo).
        :type interval: sympy Interval
        :param params: values for the kinetic parameters. All set to 1 by default.
        :type params: dictionary, key = kinetic params, value = numeric
        :rtype: sympy SparseMatrix.

        References:

        Feliu, E., & Wiuf, C. (2015). Finding the positive feedback loops underlying multi-stationarity. BMC systems biology, 9(1), 1

        """

        im = sp.SparseMatrix(self.n_species, self.n_reactions, 0)

        # variable names
        if not var:
            var = "g_"

        for j in range(self.n_reactions):
            v = self.rates[j]
            for i in range(self.n_species):
                s = self.species[i]

                # set kinetic parameters to 1
                if not params:
                    params = {}
                sbs = [f for f in v.free_symbols if str(f) not in self.species]
                for f in sbs:
                    if f not in params:
                        params[f] = 1
                v = v.subs(params)

                if check:
                    # if rate depends only on one variable, check for monotonicity
                    if len(v.atoms(sp.Symbol)) <= 1:
                        if not interval:
                            interval = sp.Interval.Lopen(0, sp.oo)
                        if not sp.is_monotonic(v, interval):
                            raise ValueError("Rate for reaction {} is not monotonic.".format(j))
                    else:
                        warnings.warn("Check for monotonicity not implemented for multivariate expressions.")

                deriv = sp.ratsimp(sp.diff(v, sp.Symbol(s)))

                # if no state is provided, use (1, 1, ..., 1)
                # and set kinetic parameters to 1
                if not state:
                    state = {}
                state = dict((f, state[f]) for f in state)
                for f in deriv.free_symbols:
                    if f not in state:
                        state[f] = 1

                # evaluate derivative in state
                deriv = deriv.subs(state)

                if deriv == sp.zoo:
                    raise ValueError("Derivative undefined for {} = {}.".format(tuple(state.keys()), tuple(state.values())))

                if deriv != 0:
                    if deriv > 0:
                        im[i, j] = sp.Symbol(var + str(i + 1) + "_" + str(j + 1))
                    else:
                        im[i, j] = -sp.Symbol(var + str(i + 1) + "_" + str(j + 1))
        return im


    ### Other attributes ###

    @property
    def deficiency(self):
        """Deficiency of the chemical reaction network,
        calculated as number of complexes, minus number of linkage classes,
        minus the rank of the stoichiometric matrix.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A <-> B + C", "2B -> C", "C -> D + E", "D + E <-> 2B"])
        >>> net.deficiency
        0

        """
        return sp.Matrix(self.incidence_matrix).rank() - sp.Matrix(self.stoich_matrix).rank()


    def format_deficiency(self):
        """Return a string 'deficiency delta = numer of complexes -
        number of linkage classes - rank of the stoichiometric matrix'

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A <-> B + C", "2B -> C", "C -> D + E", "D + E <-> 2B"])
        >>> net.format_deficiency()
        'deficiency 0 = 5 - 2 - 3'

        """
        return "deficiency {} = {} - {} - {}".format(self.deficiency, \
                                                     self.n_complexes, \
                                                     self.n_linkage_classes, \
                                                     self.stoich_matrix.rank())


    @property
    def is_ma(self):
        """Return True if the network has mass-action kinetics."""
        return all(sp.Symbol(s) not in k.cancel().atoms() for s in self.species for k in self.kinetic_params)


    @property
    def stoich_space_dim(self):
        """Return the dimension of the stoichiometric subspace,
        i.e. the rank of the stoichiometric matrix."""
        return sp.Matrix(self.stoich_matrix).rank()


    def _cons_laws(self):
        """Return a base of the nullspace in row echelon form.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A <-> B + C", "2B -> C", "C -> D + E", "D + E <-> 2B"])
        >>> net._cons_laws()
        [Matrix([
        [  1],
        [1/3],
        [2/3],
        [  0],
        [2/3]]), Matrix([
        [ 0],
        [ 0],
        [ 0],
        [ 1],
        [-1]])]
        """
        return [sp.Matrix(row) for row in sp.Matrix([list(c) for c in self.stoich_matrix.T.nullspace()]).rref()[0].tolist()]


    @property
    def cons_laws(self):
        """Tuple of generators of the network conservation laws.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A <-> B + C", "2B -> C", "C -> D + E", "D + E <-> 2B"])
        >>> net.cons_laws
        (A + B/3 + 2*C/3 + 2*E/3, D - E)
        """
        return tuple((v.T * sp.Matrix(list(map(sp.Symbol, self.species))))[0, 0] for v in self._cons_laws())


    @property
    def p_invariants(self):
        """Return a list of generators of the Petri Net P-invariants.
        Requires pycddlib."""
        return _pos_generators(self.stoich_matrix.T)


    def format_p_invariants(self):
        """Return a list of expressions in the species,
        each representing a generator of the Petri Net P-invariants.
        Requires pycddlib."""
        pinv = self.p_invariants
        if len(pinv) > 0:
            return list(pinv * sp.Matrix(list(map(sp.Symbol, self.species))))
        else:
            return []


    @property
    def t_invariants(self):
        """Return a matrix with rows the generators of the Petri Net T-invariants.
        Requires pycddlib."""
        return _pos_generators(self.stoich_matrix)


    def format_t_invariants(self):
        """Return a list of expressions in the reaction ids,
        each representing a generator of the Petri Net T-invariants.
        Requires pycddlib."""
        tinv = self.t_invariants
        if len(tinv) > 0:
            return list(tinv * sp.Matrix(list(map(sp.Symbol, self.reactionids))))
        else:
            return []


    @property
    def elem_modes(self):
        """Return the list of elementary modes of the network,
        as lists of length the number of reactions.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A -> B", "2B + 2C <-> 2D", "D -> A + C"])
        >>> net.elem_modes
        [[1, 1, 0, 1], [0, 1, 1, 0]]

        Requires pycddlib.

        """
        return self.t_invariants.tolist()


    def format_elem_modes(self):
        """Return a list of expressions in the reaction ids,
        each representing an elementary mode.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A -> B", "B + C <-> D", "D -> A + C"])
        >>> net.format_elem_modes()
        [[r0 + r1 + r2], [r1 + r1_rev]]

        Requires pycddlib.
        """
        return self.format_t_invariants()


    def is_ss_flux(self, v):
        """Check if v is a steady state flux."""
        if len(v) != self.n_reactions:
            raise ValueError("Invalid flux vector: wrong length.")
        return self.stoich_matrix * sp.Matrix(v) == sp.zeros(self.n_species, 1)


    def is_cyclic_ss_flux(self, v):
        """Check if v is a cyclic steady state flux."""
        if len(v) != self.n_reactions:
            raise ValueError("Invalid flux vector: wrong length.")
        return self.incidence_matrix * sp.Matrix(v) == sp.zeros(self.n_complexes, 1)


    def tree_constant(self, index, gma = None):
        """Return the expressions of Prop. 3 in [1], associated to *index*.
        The term ''tree constant'' is introduced in [2].

        If a dictionary gma is provided, the network is interpreted as a generalised mass action network
        with the map between complexes and kinetic complexes provided in gma,
        and the tree contant is calculated on the resulting kinetic matrix.

        :param gma: dictionary representing the map from complexes to kinetic complexes.
        :type gma: dictionary, key = string, value = string

        References:

        [1] Craciun, G. et al. (2009), Toric dynamical systems. Journal of Symbolic Computation 44.11: 1551-1565.

        [2] Johnston, M. D. (2014). Translated chemical reaction networks. Bulletin of mathematical biology, 76(5), 1081-1116.

        """
        # find the linkage classes
        slcs, lcs = self.strong_terminal_conn_components()
        l = lcs[index]
        if l not in slcs:
            warnings.warn("Complex {} not terminal.".format(self.complexes[index]))
            return None
        inds = [i for i in range(self.n_complexes) if lcs[i] == l and i != index]

        if gma:
            kparams = [r.rate / parse_complex(gma[str(r.reactant)]).ma() for r in self.reactions]
            kinetic_matrix = -self.incidence_matrix.multiply(sdiag(kparams).multiply(negative(self.incidence_matrix).T))
        else:
            kinetic_matrix = self.kinetic_matrix

        return (-1)**len(inds)*kinetic_matrix.extract(inds, inds).det()


    def tree_constants(self, gma = None):
        """Return the expressions of Prop. 3 in [1].
        The term ''tree constant'' is introduced in [2].

        If a dictionary gma is provided, the network is interpreted as a generalised mass action network
        with the map between complexes and kinetic complexes provided in gma,
        and the tree contants are calculated on the resulting kinetic matrix.

        :param gma: dictionary representing the map from complexes to kinetic complexes.
        :type gma: dictionary, key = string, value = string

        References:

        [1] Craciun, G. et al. (2009), Toric dynamical systems. Journal of Symbolic Computation 44.11: 1551-1565.

        [2] Johnston, M. D. (2014). Translated chemical reaction networks. Bulletin of mathematical biology, 76(5), 1081-1116.

        """
        if gma:
            kparams = [r.rate / parse_complex(gma[str(r.reactant)]).ma() for r in self.reactions]
            kinetic_matrix = -self.incidence_matrix.multiply(sdiag(kparams).multiply(negative(self.incidence_matrix).T))
        else:
            kinetic_matrix = self.kinetic_matrix
        # find the linkage classes
        slcs, lcs = self.strong_terminal_conn_components()
        tree_consts = []
        for index in range(self.n_complexes):
            l = lcs[index]
            if l in slcs:
                inds = [i for i in range(self.n_complexes) if lcs[i] == l and i != index]
                tree_consts.append((-1)**len(inds)*kinetic_matrix.extract(inds, inds).det())
            else:
                tree_consts.append(None)
        return tree_consts


    ### System of ODEs ###

    def equations(self):
        """Return a matrix of expressions for the derivatives of the species concentrations."""
        return self.stoich_matrix * self.rates


    def odes(self):
        """Return a list of differential equations describing the
        evolution of the species concentrations in time."""
        t = sp.Symbol('t', real = True)
        x = [sp.Function(s) for s in self.species]
        derivs = self.equations()
        return [sp.Eq(sp.Derivative(x[j](t), t),
                                    derivs[j].subs([(sp.Symbol(self.species[i]), x[i](t)) for i in range(self.n_species)]))
                for j in range(self.n_species)]


    def format_equations(self):
        """Return strings representing the differential equations
        describing the evolution of the species concentrations."""
        eqs = self.equations()
        return ["d{}/dt = {}".format(self.species[s], eqs[s, 0]) for s in range(self.n_species)]


    def groebner(self):
        """Return a groebner basis for the steady state ideal."""
        return sp.groebner([sp.ratsimp(e).as_numer_denom()[0] for e in self.equations()],
                           *[sp.Symbol(s) for s in self.species])


    def _constant_species(self):
        """Indices of species with constant concentration."""
        # This is more general than
        # [i for i in range(self.n_species) if all([j == 0 for j in self.stoich_matrix[i,:]])]
        eqs = self.equations()
        return [i for i in range(self.n_species) if eqs[i] == 0]


    @property
    def constant_species(self):
        """List of species with constant concentration."""
        return [self.species[s] for s in self._constant_species()]


    def derivative(self, cplx):
        """Return an expression for the derivative of a complex.

        cplx is a string representing the complex.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A <-> B + C", "2B -> C", "C -> D + E", "D + E <-> 2B"])
        >>> net.derivative("3A + E")
        -3*A*k_r0 + B**2*k_r3_rev + 3*B*C*k_r0_rev + C*k_r2 - D*E*k_r3
        >>> net.derivative("3A + B + 2C + 2E")
        0

        """
        c = parse_complex(cplx)
        der = 0
        # sum the derivative of each species, multiplied by the stoichiometry
        for s in c:
            if str(s) not in self.species:
                raise ValueError("Invalid complex: {} not a valid species.".format(s))
            else:
                der = der + c[s] * (self.stoich_matrix[self.species.index(str(s)),:] * self.rates)[0,0]
        return der.cancel()


    def is_constant(self, cplx):
        """Return True if the derivative of cplx is zero."""
        return self.derivative(cplx) == 0


    def has_linear_equation(self, species):
        """Check if the equation ds/dt = 0 is linear in s."""
        expr = sp.ratsimp((self.stoich_matrix[self.species.index(species), :] * self.rates)[0]).as_numer_denom()[0]
        return sp.degree(expr, sp.Symbol(species)) == 1


    def is_dyn_eq(self, net):
        """Check if two networks have the same dynamics.
        The networks must have the same set of species."""
        if set(self.species) != set(net.species):
            return False
        eqs = self.equations()
        othereqs = net.equations()
        return all((eqs[i] - othereqs[net.species.index(self.species[i])]).cancel() == 0 for i in range(self.n_species))


    ### Sources, sinks, intermediates ###

    def _simple_complexes(self):
        """Return the complexes that coincide with species."""
        return [c for c in range(self.n_complexes) \
                if sum(self.complex_matrix[:,c]) == 1 and sum(self.complex_matrix[list(self.complex_matrix[:,c]).index(1),:]) == 1]


    @property
    def simple_complexes(self):
        """Return the complexes that coincide with species."""
        return [self.complexes[c] for c in self._simple_complexes()]


    def _stoich_1_species(self):
        """Return the indices of stoichiometry 1 species."""
        return [s for s in range(self.n_species) if all([x <= 1 for x in self.complex_matrix[s,:]])]


    @property
    def stoich_1_species(self):
        """Stoichiometry 1 species."""
        return [self.species[s] for s in self._stoich_1_species()]


    def _source_complexes(self):
        """Indices of complexes that never appear as products."""
        return [c for c in range(self.n_complexes) \
                if all([i <= 0 for i in self.incidence_matrix[c,:]])]


    @property
    def source_complexes(self):
        """Complexes that never appear as products."""
        return [self.complexes[c] for c in self._source_complexes()]


    def is_source_species(self, s):
        """Check if s is a source species."""
        return all([self.complex_matrix.row(self.species.index(s)).multiply(self.incidence_matrix)[0, i] <= 0 \
                    for i in range(self.n_reactions)])


    def _source_species(self):
        """Indices of species that are never produced."""
        sm = self.stoich_matrix
        return [s for s in range(self.n_species) \
                if all([sm[s, i] <= 0 for i in range(self.n_reactions)])]


    @property
    def source_species(self):
        """Species that are never produced."""
        return [self.species[s] for s in self._source_species()]


    def _sink_complexes(self):
        """Indices of complexes that never appear as reactants."""
        return [c for c in range(self.n_complexes) \
                if all([i >= 0 for i in self.incidence_matrix[c,:]])]


    @property
    def sink_complexes(self):
        """Complexes that never appear as reactants."""
        return [self.complexes[c] for c in self._sink_complexes()]


    def is_sink_species(self, s):
        """Check if s is a sink species."""
        return all([self.complex_matrix.row(self.species.index(s)).multiply(self.incidence_matrix)[0, i] >= 0 \
                    for i in range(self.n_reactions)])


    def _sink_species(self):
        """Indices of species that are never consumed."""
        sm = self.stoich_matrix
        return [s for s in range(self.n_species) \
                if all([sm[s, i] >= 0 for i in range(self.n_reactions)])]


    @property
    def sink_species(self):
        """Species that are never consumed."""
        return [self.species[s] for s in self._sink_species()]


    def _intermediate_complexes(self):
        """Indices of complexes that are reactants for at least one reaction
        and products for at least one reaction."""
        return [c for c in range(self.n_complexes) \
                if any([i < 0 for i in self.incidence_matrix[c,:]]) and any([i > 0 for i in self.incidence_matrix[c,:]])]


    @property
    def intermediate_complexes(self):
        """Complexes that are reactants for at least one reaction
        and products for at least one reaction."""
        return [self.complexes[c] for c in self._intermediate_complexes()]


    def simple_intermediate(self, s):
        """Check whether the species is an intermediate,
        that is produced or consumed with stoichiometric coeffient 1,
        and does not interact with other species."""
        if s not in self.species:
            raise ValueError(s + " is not a valid species.")
        source = self.source_species
        sink = self.sink_species
        i = self.species.index(s)
        return all([(self.stoich_matrix.col(r)[i] == 0) or \
                    (self.stoich_matrix.col(r)[i] == 1 and all([self.stoich_matrix.col(r)[j] <= 0 for j in range(self.n_species) if j != i])) or \
                    (self.stoich_matrix.col(r)[i] == -1 and all([self.stoich_matrix.col(r)[j] >= 0 for j in range(self.n_species) if j != i]))
                    for r in range(self.n_reactions)]) and \
               s not in sink and s not in source


    def stoich_1_intermediates(self, species_list):
        """Check whether the species in species_list are intermediates,
        with stoichiometric coefficients 0 or 1."""
        for s in species_list:
            if s not in self.species:
                raise ValueError(s + " is not a valid species.")
        source = self.source_species
        sink = self.sink_species
        return all(all(i <= 1 for i in self.complex_matrix[self.species.index(s),:]) and
                   s not in sink and s not in source for s in species_list)


    def _intermediate_species(self):
        """Indices of species that are not sink or souce species."""
        source = self._source_species()
        sink = self._sink_species()
        return [s for s in range(self.n_species) if s not in source and s not in sink]


    def is_intermediate_species(self, s):
        """Check if s is an intermediate species."""
        return (not self.is_source_species(s)) and (not self.is_sink_species(s))


    @property
    def intermediate_species(self):
        """Species that are not sink or souce species."""
        return [self.species[s] for s in self._intermediate_species()]


    def _intermediate_stoich_1_species(self):
        """Indices of intermediate species that never appear with stoichiometry greater than one."""
        intermediates = self._intermediate_species()
        return [s for s in self._stoich_1_species() if s in intermediates]


    @property
    def intermediate_stoich_1_species(self):
        """Intermediate species that never appear with stoichiometry greater than one."""
        return [self.species[s] for s in self._intermediate_stoich_1_species()]


    def _intermediate_simple_complexes(self):
        """Indices of simple intermediate complexes."""
        intermediates = self._intermediate_complexes()
        simple = self._simple_complexes()
        return [c for c in intermediates if c in simple]


    @property
    def intermediate_simple_complexes(self):
        """Simple intermediate complexes."""
        return [self.complexes[c] for c in self._intermediate_simple_complexes()]


    ### Reduction ###

    def remove(self, rapid_eq = None, qss = None, cons_law = None, minimal = False, remove_const = False, \
               merge_reacts = False, adjust = False, debug = False, network_file = None):
        """Remove intermediates either by rapid-equilibrium or
        by quasi-steady state approximation, with an optional conservation law.

        :param rapid_eq: list of pairs species-complexes. The species are to be replaced by complexes using rapid equilibrium approximation.
        :type rapid_eq: list of (string (species), string (complex)).
        :param qss: list of species to remove via qssa.
        :type qss: list of strings
        :param cons_law: remove a species using a conservation. All other species involved in
                         in the conservation are eliminated first, via quasi-steady state approximation,
                         if not listed in *rapid_eq*.
        :type cons_law: (string, ConsLaw)
        :param minimal: find network of minimal structure when applying qss.
        :type minimal: boolean
        :param remove_const: remove any species left constant after the reduction.
        :type remove_const: boolean
        :param merge_reacts: merge reactions with the same reactant and product.
        :type merge_reacts: boolean
        :param adjust: change the rates so that they reflect the reactant species.
        :type adjust: boolean
        :param network_file: save the reduction steps to file with the given path.
        :type network_file: string

        Update remove_species.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["s + e (k_1)<->(k1) es", "es ->(k2) e + p", "i + e (k_3)<->(k3) ei", "i + es (k_3)<->(k3) esi", "s + ei (k_1)<->(k1) esi"])
        >>> cl = ('e', ConsLaw('e + ei + es + esi', 'et'))
        >>> net.remove(rapid_eq = [('ei', 'e + i'), ('esi', 'e + s + i'), ('es', 's + e')], cons_law = cl)
        >>> net.reactions
        (r1: s ->(et*k1*k2*k_3/(i*k1*k3*s + i*k3*k_1 + k1*k_3*s + k_1*k_3)) p,)

        """

        if cons_law:
            add_species = [str(c) for c in cons_law[1].species.keys() if c != sp.Symbol(cons_law[0])]
        else: add_species = []

        if rapid_eq != None: rapid_eq_species = [pair[0] for pair in rapid_eq]
        else: rapid_eq_species, rapid_eq = [], []
        if qss != None: qss_species = qss + [a for a in add_species if a not in rapid_eq_species and a not in qss]
        else: qss_species = [a for a in add_species if a not in rapid_eq_species]

        if network_file: self.save_to_file(network_file, overwrite = 'w', log = "Original network")

        # rapid-equilibrium
        for pair in rapid_eq:
            self._rapid_eq(pair, debug)
            if network_file: self.save_to_file(network_file, log = "Rapid eq. on {}, {}".format(pair[0], pair[1]))
        if debug: print("removed_species after rapid_eq: {}".format(self.removed_species))

        # qss
        self._qss(qss_species, minimal = minimal, network_file = network_file, debug = debug)
        if debug: print("removed_species after QSS: {}".format(self.removed_species))

        # use conservation law to rewrite cons_law[0] in terms of remaining variables
        if cons_law:
            self.remove_by_cons(cons_law[0], cons_law[1], debug)

        # optional changes
        if remove_const:
            self.remove_all_constants()
            if network_file: self.save_to_file(network_file, log = "Removed constants")
        if merge_reacts:
            self.merge_reactions()
            if network_file: self.save_to_file(network_file, log = "Merged reactions")
        if adjust: self._fix_ma()

        if network_file: self.save_to_file(network_file, rs = True, log = "Final network")


    def remove_by_cons(self, species, cons_law, debug = False):
        """Remove a species using a conservation law.
        First replace removed_species in the conservation law with their expression.
        Then use the conservation expression to write the species
        concentration as function of the remaining species.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["E + S (k_1)<->(k1) C", "C ->(k2) E + P"])
        >>> net.qss("C")
        >>> net.reactions
        (r0_r1: E + S ->(k1*k2/(k2 + k_1)) E + P,)
        >>> net.removed_species
        (('C', E*S*k1/(k2 + k_1)),)
        >>> cl = ConsLaw("E + C", "etot")
        >>> net.remove_by_cons("E", cl)
        >>> net.reactions
        (r0_r1: S ->(etot*k1*k2/(S*k1 + k2 + k_1)) P,)

        References:

        Tonello et al. (2016), On the elimination of intermediate species in chemical reaction networks.

        """
        conservation = cons_law.expression

        if debug:
            print("Removed species: {}".format(self.removed_species))
            print("Conservation: {}".format(conservation))

        for variable, expr in self.removed_species:
            if debug:
                print("Replacing {} with {}".format(variable, expr))
            conservation = conservation.subs(variable, expr)
            if debug:
                print("Found {}".format(conservation))
                print

        # The next is quicker, but not always applicable
        #conservation = (conservation / sp.Symbol(species)).cancel()
        #exp = cons_law.constant / conservation
        exp = sp.solve(conservation - cons_law.constant, sp.Symbol(species))[0]

        # remove species
        self.remove_constant(species, expr = exp)

        if debug: print("Remove by Conservation: added to removed_species {}".format(self.removed_species))


    def qss(self, intermediate = None, cons_law = None, minimal = False, remove_const = False, \
            merge_reacts = False, adjust = False, debug = False, network_file = None):
        """Remove an intermediate species via quasi-steady state approximation.

        Keyword arguments:

        :param intermediate: species to remove via qssa.
        :type intermediate: string
        :param cons_law: remove a species using a conservation. All other species involved in
                         in the conservation are eliminated first, via quasi-steady state approximation.
        :type cons_law: (string, ConsLaw)
        :param minimal: find network of minimal structure when applying qss.
        :type minimal: boolean
        :param remove_const: remove any species left constant after the reduction.
        :type remove_const: boolean
        :param merge_reacts: merge reactions with the same reactant and product.
        :type merge_reacts: boolean
        :param adjust: change the rates so that they reflect the reactant species.
        :type adjust: boolean
        :param network_file: save the reduction steps to file with the given path.
        :type network_file: string

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["E + S (k_1)<->(k1) C", "C ->(k2) E + P"])
        >>> net.qss("C")
        >>> net.reactions
        (r0_r1: E + S ->(k1*k2/(k2 + k_1)) E + P,)
        >>> net.removed_species
        (('C', E*S*k1/(k2 + k_1)),)
        >>> net = from_react_strings(["E + S (k_1)<->(k1) C", "C ->(k2) E + P"])
        >>> net.qss("C", cons_law = ("E", ConsLaw("E + C", "etot")))
        >>> net.reactions
        (r0_r1: S ->(etot*k1*k2/(S*k1 + k2 + k_1)) P,)


        Update removed_species.

        References:

        Tonello et al. (2016), On the elimination of intermediate species in chemical reaction networks.

        Madelaine et al. (2016), Normalizing Chemical Reaction Networks by Confluent Structural Simplification (CMSB).

        """
        return self.remove(qss = ([intermediate] if intermediate else None), cons_law = cons_law,
                           adjust = adjust, minimal = minimal, debug = debug, \
                           network_file = network_file, remove_const = remove_const, merge_reacts = merge_reacts)


    def _qss(self, intermediates, minimal = False, error_if_missing = True, network_file = None, debug = False):
        """Eliminate the intermediates via quasi-steady state approximation."""

        if debug: print("Intermediates to remove: {}".format(intermediates))
        if minimal:
            fvs = dict((self.reactions[r].reactionid, \
                        [0 if i != r else 1 for i in range(self.n_reactions)]) \
                        for r in range(self.n_reactions)) # flux vectors
            nr = self.n_reactions # flux vector length

        reactions = self.reactions

        for intermediate in intermediates:
            # if intermediate not in species, raise error or warning
            if intermediate not in self.species:
                if error_if_missing:
                    raise ValueError("Species {} not in network.".format(intermediate))
                else:
                    warnings.warn("Species {} not in network.".format(intermediate))
                    return

            if debug: print("Removing: {}".format(intermediate))

            # Error if species has stoichiometry greater than 1
            if any(s > 1 for s in self.complex_matrix[self.species.index(intermediate), :]):
                raise ValueError("Species {} appears with stoichiometry greater than 1.".format(intermediate))

            # Error if species has a non linear derivative
            hasLinearDyn = self.has_linear_equation(intermediate)
            if not hasLinearDyn:
                raise ValueError("Species {} has nonlinear kinetics.".format(intermediate))

            gens = [] # generators
            newreactions = []
            reactReactions = []
            prodReactions = []

            # case of intermediate in both reactant and product:
            # move intermediate to rate
            for r in reactions: r.remove_react_prod(intermediate)

            # separate reactions into reactions producing intermediate
            # and reactions consuming the intermediate
            for r in range(len(reactions)):
                if intermediate in reactions[r].product: reactReactions.append(reactions[r])
                else:
                    if intermediate in reactions[r].reactant: prodReactions.append(reactions[r])
                    else:
                        newreactions.append(reactions[r])
                        # add flux vector to generators
                        if minimal:
                            gens.append(fvs[reactions[r].reactionid])

            y = sp.Symbol(intermediate)
            expr = y
            # if there are no reactant or product reactions, species is not an intermediate
            if len(prodReactions) == 0 or len(reactReactions) == 0:
                raise ValueError("Species {} is not a valid intermediate.".format(intermediate))

            # find expression for intermediate
            # (does not work in more general cases)
            numer = sum([reaction.rate * reaction.product[intermediate] for reaction in reactReactions])
            denom = sum([reaction.reactant[intermediate] * reaction.rate for reaction in prodReactions])
            expr = (y * numer / denom).cancel()

            def combine(r1, r2):
                newid = r1.reactionid + '_' + r2.reactionid
                c = Complex(dict((k, d) for k, d in r1.product.items() if k != intermediate))
                cprime = Complex(dict((k, d) for k, d in r2.reactant.items() if k != intermediate))
                cprimeprime = c | cprime
                newreactant = Complex(r1.reactant + cprimeprime - c)
                newproduct = Complex(r2.product + cprimeprime - cprime)
                return Reaction(newid, newreactant, newproduct, (r1.rate * r2.rate / denom).cancel())

            # create reactions by combining reactions
            # producing and consuming the intermediate
            for r1 in reactReactions:
                for r2 in prodReactions:
                    coeffs = None
                    if minimal:
                        # flux vector of the combination
                        v = [fvs[r1.reactionid][i] + fvs[r2.reactionid][i] for i in range(nr)]
                        # check if v is generator
                        coeffs = _pos_dependent(gens, v)

                    if coeffs:
                        if debug:
                            print("{}, {} not a generator".format(r1.reactionid, r2.reactionid))
                        for i in range(len(gens)):
                            if coeffs[i] != 0:
                                newreactions[i]._rate = (newreactions[i].rate + \
                                                        r1.rate * r2.rate / denom * \
                                                        sp.nsimplify(coeffs[i], rational = True)).cancel()
                    else:
                        newreactions.append(combine(r1, r2))
                        if minimal:
                            fvs[newreactions[-1].reactionid] = v
                            gens.append(v)

            # replace all instances of intermediate with expression
            for r in range(len(newreactions)):
                newreactions[r]._rate = (newreactions[r].rate).subs(y, expr).factor()

            self._removed_species.append((intermediate, expr))
            reactions = list(newreactions)
            if network_file:
                self.save_to_file(network_file,
                                  reactions = [r for r in reactions if r.reactant != r.product],
                                  log = "Qss on " + intermediate)

        # remove cycles
        reactions = [r for r in reactions if r.reactant != r.product]

        if debug:
            print("New reactions")
            for r in reactions: print(r)
            print

        self.reactions = reactions


    def _qss_generalised(self, intermediate, keep_loops = False,
                         no_rates = False, error_if_missing = True, debug = False):
        """Implement the quasi-steady state approximation for stoich > 1."""

        # if intermediate not in species, raise error or warning
        if intermediate not in self.species:
            if error_if_missing:
                raise ValueError("Species {} not in network.".format(intermediate))
            else:
                warnings.warn("Species {} not in network.".format(intermediate))
                return

        # if species is constant, remove it from complexes
        if self.is_constant(intermediate):
            self.remove_constant(intermediate)
            return

        # Notify if species has a non linear derivative
        hasLinearDyn = self.has_linear_equation(intermediate)
        if not hasLinearDyn:
            raise ValueError("Species {} has nonlinear kinetics.".format(intermediate))

        newreactions = []
        reactReactions = []
        prodReactions = []

        # case of intermediate in both reactant and product:
        # move intermediate to rate
        reactions = self.reactions
        for r in reactions: r.remove_react_prod(intermediate)

        # separate reactions to those producing intermediate
        # and those using the intermediate
        for r in range(len(reactions)):
            if intermediate in reactions[r].product:
                reactReactions.append(reactions[r])
            else:
                if intermediate in reactions[r].reactant:
                    prodReactions.append(reactions[r])
                else:
                    newreactions.append(reactions[r])

        y = sp.Symbol(intermediate)
        expr = y
        # if there are no reactant or product reactions, species is not an intermediate
        if len(prodReactions) == 0 or len(reactReactions) == 0:
            raise ValueError("Species is not a valid intermediate.")
        else:
            if not no_rates:
                if hasLinearDyn:
                    denom = sum([reaction.reactant[intermediate] * reaction.rate for reaction in prodReactions])
                    expr = sp.solve(sp.ratsimp((self.stoich_matrix[self.species.index(intermediate), :] * self.rates)[0]).as_numer_denom()[0], y)[0]
                else:
                    expr = sp.solve((self.stoich_matrix[self.species.index(intermediate), :] * self.rates)[0], y)[1]
                    denom = 1

            for r1 in reactReactions:
                for r2 in prodReactions:
                    n = r1.product[intermediate]
                    m = r2.reactant[intermediate]
                    h = sp.lcm(n, m)
                    c = Complex(dict((k, d) for k, d in r1.product.items() if k != intermediate))
                    cprime = Complex(dict((k, d) for k, d in r2.reactant.items() if k != intermediate))
                    hn = int(h / n)
                    hm = int(h / m)
                    cprimeprime = c.times(hn) | cprime.times(hm)

                    newreactant = Complex(r1.reactant.times(hn) + cprimeprime - c.times(hn))
                    newproduct = Complex(r2.product.times(hm) + cprimeprime - cprime.times(hm))
                    if keep_loops or (newreactant != newproduct):
                        # we add "_" in front if there is a number
                        # because a reaction id can not start with a number
                        newid = ("" if h / n == 1 else "_" + str(h / n)) + r1.reactionid + '_' + \
                                ("" if h / m == 1 else str(h / m)) + r2.reactionid
                        newreactions.append(Reaction(newid, \
                                                     newreactant, \
                                                     newproduct, \
                                                     sp.Symbol("k_" + newid) * newreactant.ma() if no_rates
                                                     else (sp.gcd(n, m) * r1.rate * r2.rate / denom).cancel()))

            for r in range(len(newreactions)):
                newreactions[r]._rate = (newreactions[r].rate).subs(y, expr).cancel()

        if debug:
            print("New reactions")
            for r in newreactions: print(r)
            print

        self.reactions = newreactions
        self._removed_species.append((intermediate, expr))


    def _fix_ma(self):
        """Adjust the reactant and products so that they match the
        monomials in the numerators of the reaction rates."""
        reactions = list(itertools.chain(*map(lambda r: _split_reaction_monom(r, self.species), self.reactions)))
        for r in reactions: r._fix_ma(self.species)
        self.reactions = reactions


    def _remove_react_prod(self):
        """For each reaction, remove the species that
        are in common between reactant and product."""
        for r in self.reactions: r.remove_react_prod()
        self.reactions = self.reactions


    def _same_denom(self):
        self.reactions = _same_denom(self.reactions)


    def _fix_denom(self):
        """Remove reactant and product species if their
        concentrations divides the denominator of the rate."""
        reactions = list(itertools.chain(*map(lambda r: _split_reaction_monom(r, self.species), self.reactions)))
        for r in reactions: r._fix_denom(self.species)
        self.reactions = reactions


    def rapid_eq_with_pool(self, remove, keep, pool_name = None, cons_law = None, debug = False):
        """Apply rapid equilibrium to reactions between remove and keep.
        remove and keep are merged in a pool."""
        # standardise complexes
        remove = parse_complex(remove)
        keep = parse_complex(keep)

        # Check that complexes exist
        if not remove in self.complexes:
            raise ValueError(remove + " is not a valid complex.")
        if not keep in self.complexes:
            raise ValueError(keep + " is not a valid complex.")

        # Index of complex to be removed and of remaining complex
        removeComplex = self.complexes.index(remove)
        remove = str(remove)
        removeSpecies = self.species.index(remove)

        keepComplex = self.complexes.index(keep)
        keepSpecies = [s for s in range(self.n_species) if self.complex_matrix[s, keepComplex] != 0]

        # Check that forward and backward reactions exist
        # and find index
        forward = None
        backward = None
        for r in range(self.n_reactions):
            if self.incidence_matrix.col(r)[removeComplex] == -1 and self.incidence_matrix.col(r)[keepComplex] == 1:
                forward = r
            if self.incidence_matrix.col(r)[removeComplex] == 1 and self.incidence_matrix.col(r)[keepComplex] == -1:
                backward = r

        if forward == None:
            raise ValueError("No reaction from {} to {}.".format(remove, keep))

        if backward == None:
            raise ValueError("No reaction from {} to {}.".format(keep, remove))

        exp = sp.solve(self.rates[forward] - self.rates[backward], sp.Symbol(remove))[0]
        self.logger.info('Found ' + remove + ' = ' + str(exp))

        if pool_name == None:
            pool_name = "pool_" + remove + "_" + keep
        pool_complex = parse_complex(pool_name)

        # find pool is keep + remove
        # find remove and keep as functions of pool
        keep_pool = sp.Symbol(pool_name) / (1 + exp / sp.Symbol(str(keep)))
        keep_factor = (1 / (1 + exp / sp.Symbol(str(keep)))).cancel()
        remove_pool = keep_pool * exp / sp.Symbol(str(keep))
        remove_factor = (keep_factor * exp / sp.Symbol(str(keep))).cancel()

        reactions = []
        for reaction in self.reactions:
            reactant = reaction.reactant
            product = reaction.product
            rate = reaction.rate
            if reaction.product == keep or str(reaction.product) == str(remove):
                product = pool_complex
            else:
                product = reaction.product
            if reaction.reactant == keep:
                reactant = pool_complex
                rate = (keep_factor * rate * pool_complex.ma() / sp.Symbol(str(keep))).cancel()
            if str(reaction.reactant) == remove:
                reactant = pool_complex
                rate = (remove_factor * rate * pool_complex.ma() / sp.Symbol(remove)).cancel()
            if reactant != product:
                reactions.append(Reaction(reaction.reactionid, \
                                          reactant, \
                                          product, \
                                          rate.cancel()))

        self.reactions = reactions
        self._removed_species.append((remove, exp))

        if debug:
            for r in self.reactions: print(r)

        if cons_law:
            self.remove_by_cons(cons_law[0], cons_law[1], debug)


    def rapid_eq(self, species, cmplx, cons_law = None, debug = False, network_file = None):
        """Apply the rapid equilibrium approximation to the reaction from species to cmplx,
        replacing the species with the complex cmplx.

        Keyword arguments:

        :param species: species to eliminate.
        :type species: string
        :param complex: complex supposed at rapid equilibrium with species.
        :type complex: string
        :param cons_law: remove a species using a conservation. All other species involved in
                         in the conservation are eliminated first, with the exception of *species*,
                         via quasi-steady state approximation.
        :type cons_law: (string, ConsLaw)
        :param network_file: save the reduction steps to file with the given path.
        :type network_file: string

        Update remove_species.

        :Example:

        >>> net = from_react_strings(["E + S (k_1)<->(k1) C", "C ->(k2) E + P"])
        >>> net.rapid_eq(("C", "E + S"))
        >>> net.reactions
        (r1: E + S ->(k1*k2/k_1) E + P,)
        >>> net.removed_species
        (('C', E*S*k1/k_1),)
        >>> net = from_react_strings(["E + S (k_1)<->(k1) C", "C ->(k2) E + P"])
        >>> net.rapid_eq(("C", "E + S"), cons_law = ("E", ConsLaw("E + C", "etot")))
        >>> net.reactions
        (r1: S ->(etot*k1*k2/(S*k1 + k_1)) P,)

        References:

        Tonello et al. (2016), On the elimination of intermediate species in chemical reaction networks.

        """
        return self.remove(rapid_eq = [(species, cmplx)], cons_law = cons_law, debug = debug, network_file = network_file)


    def _rapid_eq(self, recompls, debug):
        """The reactions between the complexes in recompls
        are assumed at rapid equilibrium.
        The first complex is removed and the second
        replaces it in the network."""

        remove, keep = recompls
        removecomplex, keepcomplex = map(parse_complex, recompls)

        # check that remove is a valid species
        if remove not in self.species:
            raise ValueError("{} not a valid species.".format(remove))
        # check that keep is a valid complex
        if keepcomplex not in self.complexes:
            raise ValueError("{} not a valid complex.".format(keep))

        reactions = self.reactions

        # Check that forward and backward reactions exist
        # and find index
        forward = None
        backward = None
        for r in range(self.n_reactions):
            if reactions[r].reactant == keepcomplex and reactions[r].product == removecomplex:
                forward = r
            if reactions[r].reactant == removecomplex and reactions[r].product == keepcomplex:
                backward = r

        if forward == None:
            raise ValueError("No reaction from {} to {}.".format(remove, keep))

        if backward == None:
            raise ValueError("No reaction from {} to {}.".format(keep, remove))

        eq = sp.ratsimp(reactions[forward].rate - reactions[backward].rate).as_numer_denom()[0]
        if sp.degree(eq, sp.Symbol(remove)) != 1:
            raise ValueError("Unable to remove {} via rapid equilibrium.".format(remove))
        exp = sp.solve(eq, sp.Symbol(remove))
        if len(exp) < 1 or exp[0] == 0:
            raise ValueError("Unable to remove {} via rapid equilibrium.".format(remove))
        exp = exp[0]
        self.logger.info('Found ' + remove + ' = ' + str(exp))

        newreactions = []
        for reaction in reactions:
            newreactant = Complex(dict(reaction.reactant))
            newproduct = Complex(dict(reaction.product))
            newrate = reaction.rate
            if remove in reaction.reactant:
                newreactant = Complex(newreactant - \
                                      removecomplex.times(reaction.reactant[remove]) + \
                                      keepcomplex.times(reaction.reactant[remove]))
            if remove in reaction.product:
                newproduct = Complex(newproduct - \
                                     removecomplex.times(reaction.product[remove]) + \
                                     keepcomplex.times(reaction.product[remove]))
            if newreactant != newproduct:
                newreactions.append(Reaction(reaction.reactionid, \
                                             newreactant, \
                                             newproduct, \
                                             newrate))

        for r in range(len(newreactions)):
            newreactions[r]._rate = (newreactions[r].rate).subs(sp.Symbol(remove), exp)

        self.reactions = newreactions

        self._removed_species.append((remove, exp))
        if debug: print("Rapid Eq: added to removed_species {}".format(self.removed_species))


    def remove_constant(self, const_species, expr = None, debug = False):
        """Remove the species from the network by setting it to constant.
        Give a warning if the derivative of the species concentration
        is not constant in the original network.

        :param expr: optional expression that replaces *const_species* in rates.
        :type expr: sympy expression

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["E + S (k_1)<->(k1) E + P"])
        >>> net.remove_constant("E", "etot")
        >>> net.reactions
        (r0: S ->(etot*k1) P, r0_rev: P ->(etot*k_1) S)

        """

        # check that const_species is a valid species
        # either const_species is in set of species, and has null derivative
        # or is used in kinetics
        if const_species not in self._species + [str(v) for v in list(set(list(itertools.chain(*[k.atoms() for k in self.kinetic_params]))))]:
            warnings.warn("Species {} not valid.".format(const_species))
        else:
            if const_species in self._species and self.derivative(const_species) != 0:
                warnings.warn("Concentration of {} not constant.".format(const_species))

        self.logger.info("Removing: {}".format(', '.join(const_species)))

        # remove species from reactants and products
        # and replace species with expression in rates
        for reaction in self.reactions:
            if const_species in reaction.reactant:
                del reaction.reactant[const_species]
                reaction._rate = reaction.rate
            if const_species in reaction.product:
                del reaction.product[const_species]
            if expr != None: reaction._rate = reaction.rate.subs(sp.Symbol(const_species), expr).factor()

        reactions = [r for r in self.reactions if r.reactant != r.product]
        self.reactions = reactions
        self._removed_species = self._removed_species + [(const_species, expr if expr else sp.Symbol(const_species))]


    def remove_all_constants(self, debug = False):
        """Remove all constant species."""
        for x in self.constant_species:
            self.remove_constant(x, debug = debug)


    def merge_reactions(self):
        """Merge the reactions with same reactant and same product."""
        self.reactions = merge_reactions(self.reactions)


    ### Graphs ###

    def complex_graph_adj(self):
        """Return the adjacency matrix of the graph of complexes.

        This is a sympy Matrix of dimension number of complexes times number of complexes.
        """
        adjacency = - self.incidence_matrix.multiply(negative(self.incidence_matrix).T)
        for i in range(self.n_complexes): adjacency[i, i] = 0
        return adjacency


    def dsr_graph_adj(self, keep = None):
        """Return the adjacency matrix of the directed species-reaction graph.

        Optionally remove the variables of the influence matrix not in *keep*.

        References:

        Feliu, E., & Wiuf, C. (2015). Finding the positive feedback loops underlying multi-stationarity. BMC systems biology, 9(1), 1

        """
        A = self.stoich_matrix
        im = self.influence_matrix()
        if keep is not None:
            im = im.applyfunc(lambda x: x if x in keep else 0)
        im = im.applyfunc(lambda x: -1 if str(x)[0] == "-" else (1 if x != 0 else 0))
        adjacency = sp.zeros(im.rows, im.rows).row_join(im).col_join(A.T.row_join(sp.zeros(A.cols, A.cols)))
        return adjacency


    ### Connectivity properties ###

    def strong_conn_components(self):
        """Return the number of strongly connected components of the graph of complexes,
        and an array of labels for the strongly connected components."""
        adjacency = csr_matrix(np.array(self.complex_graph_adj().tolist()).astype(np.int))
        return connected_components(adjacency, connection = 'strong')


    def strong_terminal_conn_components(self):
        """Return the list labels for the terminal strongly connected components,
        of the graph of complexes, and the array of labels for the strongly connected components."""
        n, cs = self.strong_conn_components()
        adjacency = csr_matrix(np.array(self.complex_graph_adj().tolist()).astype(np.int))
        stcc = []
        for i in range(n):
            complexes = [j for j in range(self.n_complexes) if cs[j] == i]
            reached = [h for h in range(self.n_complexes) if any(adjacency[h, j] == 1 for j in complexes)]
            if all([cs[h] == i for h in reached]): stcc.append(i)
        return stcc, cs


    def _terminal_complexes(self):
        """Return the indices of the complexes in the
        terminal strongly connected components."""
        stcc, cs = self.strong_terminal_conn_components()
        return [c for c in range(self.n_complexes) if cs[c] in stcc]


    @property
    def terminal_complexes(self):
        """Return the list of complexes in the terminal strongly connected components."""
        return [self.complexes[c] for c in self._terminal_complexes()]


    def _non_terminal_complexes(self):
        """Return the indices of the complexes in the
        nonterminal strongly connected components."""
        stcc, cs = self.strong_terminal_conn_components()
        return [c for c in range(self.n_complexes) if cs[c] not in stcc]


    @property
    def non_terminal_complexes(self):
        """Return the list of complexes in the nonterminal strongly connected components."""
        return [self.complexes[c] for c in self._non_terminal_complexes()]


    def weak_conn_components(self):
        """Return the number of weakly connected components of the graph of complexes,
        and an array of labels for the weakly connected components."""
        return connected_components(csr_matrix(np.array(self.complex_graph_adj().tolist()).astype(np.int)))


    @property
    def n_linkage_classes(self):
        """Number of linkage classes of the graph of complexes."""
        return self.n_complexes - sp.Matrix(self.incidence_matrix).rank()


    @property
    def linkage_classes(self):
        """List of complexes grouped by linkage class.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A <-> B + C", "2B -> C", "C -> D + E", "D + E <-> 2B"])
        >>> net.linkage_classes
        [[A, B + C], [2B, C, D + E]]

        """
        nlc, lcs = self.weak_conn_components()
        return [[self.complexes[c] for c in range(self.n_complexes) if lcs[c] == lc] for lc in range(nlc)]


    @property
    def strong_linkage_classes(self):
        """List of complexes grouped by strong linkage class.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A <-> B + C", "2B -> C", "C -> D + E", "D + E <-> 2B", "C -> F"])
        >>> net.strong_linkage_classes
        [[A, B + C], [2B, C, D + E], [F]]

        """
        nlc, lcs = self.strong_conn_components()
        return [[self.complexes[c] for c in range(self.n_complexes) if lcs[c] == lc] for lc in range(nlc)]


    @property
    def is_rev(self):
        """Return True if for each reaction with reactant c and product c' in the network
        there is at least one reaction with reactant c' and product c."""
        adj = self.complex_graph_adj().applyfunc(lambda x: 1 if x else 0)
        return adj.is_symmetric(simplify = False)


    @property
    def is_weakly_rev(self):
        """Return True if the network is weakly reversible."""
        return self.strong_conn_components()[0] == self.weak_conn_components()[0]


    ### Robustness ###

    def acr_species(self, subnets = False, same_ems = False):
        """Return some of the species with absolute concentration robustness,
        using the methods for deficiency 0 and deficiency 1 networks.

        * subnets -- if True, look for a decomposition in subnetworks using split_by_ems,
                   and use the decomposition to look for robust ratios and species.
        * same_ems -- if True and subnets = True, use acr_same_ems.

        :Example:

        >>> net = from_react_strings(["A + B -> 2B", "B -> A", "2A <-> C", "A + C <-> D"])
        >>> net.acr_species()
        ['A']
        >>> net.acr_species(subnets=True)
        ['A', 'C', 'D']

        References:

        Shinar, G., Feinberg, M. (2010). Structural sources of robustness in biochemical reaction networks. Science.

        """
        acr_s = []
        acr_c = self.acr_complexes(as_vectors = True, subnets = subnets, same_ems = same_ems)
        rank = sp.Matrix(acr_c).rank()
        if acr_c:
            for i in range(self.n_species):
                v = sp.Matrix([0 if j != i else 1 for j in range(self.n_species)]).T
                if sp.Matrix(acr_c).col_join(v).rank() == rank:
                    acr_s.append(self.species[i])
        return acr_s


    def acr_complexes(self, as_vectors = False, subnets = False, same_ems = False):
        """Return ratios of monomials that are robust,
        i.e. ratios that take the same value at each positive steady state,
        using the methods for deficiency 0 and deficiency 1 networks.

        * as_vectors -- if True, return a list of vectors of length = n_species, defining the ratios.
        * subnets -- if True, look for a decomposition in subnetworks using split_by_ems,
                     and use the decomposition to look for robust ratios.
        * same_ems -- if True and subnets = True, use acr_same_ems.

        :Example:

        >>> net = from_react_strings(["A + B -> 2B", "B -> A", "2A <-> C", "A + C <-> D"])
        >>> net.acr_complexes()
        [A]
        >>> net.acr_complexes(as_vectors=True)
        [[1, 0, 0, 0]]
        >>> net.acr_complexes(subnets=True)
        [A**2/C, A*C/D, A]

        References:

        Shinar, G., Feinberg, M. (2010). Structural sources of robustness in biochemical reaction networks. Science.

        """
        if self.is_ma:
            if subnets:
                acr_cs = self.acr_same_ems(as_vectors = True) if same_ems else []
                nets = self.subnets()
                for net in nets:
                    convert = lambda v: [v[net.species.index(self.species[i])]
                                           if self.species[i] in net.species else 0
                                           for i in range(len(self.species))]
                    acr_cs = acr_cs + list(map(convert, net.acr_complexes(as_vectors = True)))
                acr_cs = list(map(list, sp.Matrix(acr_cs).T.columnspace()))
                if as_vectors:
                    return acr_cs
                else:
                    return list(map(self._vect_to_monom, acr_cs))
            else:
                deficiency = self.deficiency
                wr = self.is_weakly_rev
                if (deficiency == 0 and wr) or deficiency == 1:
                    if deficiency == 1:
                        indgroups = [self._non_terminal_complexes()]

                    if deficiency == 0 and wr:
                        nlc, lc = self.weak_conn_components()
                        indgroups = [[i for i in range(self.n_complexes) if lc[i] == l] for l in range(nlc)]

                    m = list(map(list, sp.Matrix([list(self.complexes[i].to_vector(self.species) -
                                                       self.complexes[j].to_vector(self.species))
                                                  for inds in indgroups
                                                  for i, j in itertools.combinations(inds, 2)]).T.columnspace()))

                    if as_vectors:
                        return m
                    else:
                        return list(map(self._vect_to_monom, m))
        return []


    def _vect_to_monom(self, v):
        return sp.Mul(*(sp.Symbol(self.species[i])**v[i] for i in range(self.n_species)))


    def acr_same_ems(self, as_vectors = False):
        """Return ratios of monomials that are robust,
        i.e. ratios that take the same value at each positive steady state,
        by identifying reactions that take part in the same elementary modes,
        with the same coefficients.

        * as_vectors -- if True, return a list of vectors of length = n_species, defining the ratios.
        """
        if self.is_ma:
            tinvs = self.t_invariants.T
            if tinvs:
                # group rows
                rinds = defaultdict(list)
                for r in range(tinvs.rows):
                    rinds[tuple(tinvs.row(r))].append(r)
                cgroups = [[self.reactions[r].reactant for r in v] for v in rinds.values() if len(v) > 1]
                m = list(map(list, sp.Matrix([list(cs[i].to_vector(self.species) - cs[j].to_vector(self.species))
                                              for cs in cgroups
                                              for i, j in itertools.combinations(range(len(cs)), 2)]).T.columnspace()))
                if as_vectors:
                    return m
                else:
                    return list(map(self._vect_to_monom, m))
        return []


    ### Subnetworks ###

    def subnets(self):
        """Try to split the network into subnetworks such that
        the ranks of the subnetworks sum to the rank of the network.
        Use split_by_ems to look for candidate subnetworks.

        :rtype: list of CRN objects.
        """
        reacts, rem_reacts = self.split_by_ems()
        nets = [from_reacts(r) for r in reacts] + ([from_reacts(rem_reacts)] if rem_reacts else [])
        if sum([net.stoich_matrix.rank() for net in nets]) == self.stoich_matrix.rank():
            return nets
        else:
            return [self]


    def split_by_ems(self, same_react = False, warn = False):
        """Split reactions according to the elementary modes they take part in.

        Return a list of reactions grouped by elementary mode, and the list of reactions
        not taking part in any elementary mode.

        :param same_react: if True, do not split reactions with the same reactant.
        :type same_react: boolean
        :param warn: if True, give a warning if not all reactions take part in at least
                     one elementary mode.
        :type warn: boolean

        :rtype: (list of lists reactions, list of reactions).
        """
        subnet = {}
        tinvs = self.t_invariants

        if tinvs.rows == 0:
            return [self.reactions], []

        # case of reactions not taking part in any em
        if warn:
            if any(sum(tinvs[:, c]) == 0 for c in range(tinvs.cols)):
                warnings.warn("Network not covered by elementary modes.")

        a = sp.SparseMatrix(sp.zeros(self.n_reactions))
        for t in range(tinvs.rows):
            inds = [r for r in range(self.n_reactions) if tinvs[t, r] != 0]
            for i in inds:
                for j in inds:
                    a[i, j] = 1
        if same_react:
            for c in self.complexes:
                inds = [r for r in range(self.n_reactions) if self.reactions[r].reactant == c]
                for i in inds:
                    for j in inds:
                        a[i, j] = 1
        ncc, cc = connected_components(csr_matrix(np.array(a.tolist()).astype(np.int)))
        rcc = [[self.reactions[r] for r in range(self.n_reactions) if cc[r] == l] for l in range(ncc)]

        return [rc for rc in rcc if len(rc) > 1], [rc[0] for rc in rcc if len(rc) == 1]


    ### Print ###

    def print_laplacian(self, numeric = False, precision = 3):
        """Print the laplacian, with labels for the complexes.

        If numeric = True, the matrix is converted to a numpy array,
        and precision can be used to format the elements (default is 3).
        """
        print_matrix(self.laplacian,
                     list(map(str, self.complexes)),
                     list(map(str, self.complexes)),
                     numeric,
                     precision)


    def print_kinetic_matrix(self, numeric = False, precision = 3):
        """Print the kinetic matrix, with labels for the complexes.

        If numeric = True, the matrix is converted to a numpy array,
        and precision can be used to format the elements (default is 3).
        """
        print_matrix(self.kinetic_matrix,
                     list(map(str, self.complexes)),
                     list(map(str, self.complexes)),
                     numeric,
                     precision)


    def print_stoich_matrix(self):
        """Print the stoichiometrix matrix, with labels for species and reactions.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 (k_2)<->(k2) 2 A3", "A3 ->(k3) "])
        >>> net.print_stoich_matrix()
             r0  r1  r1_rev  r2
        A1 | -1   0       0   0 |
        A2 |  1  -1       1   0 |
        A3 |  1   2      -2  -1 |

        """
        print_matrix(self.stoich_matrix, self.species, self.reactionids)


    def print_incidence_matrix(self):
        """Print the incidence matrix, with labels for complexes and reactions."""
        print_matrix(self.incidence_matrix, list(map(str, self.complexes)), self.reactionids)


    def print_complex_matrix(self):
        """Print the complex matrix, with labels for species and complexes.

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 (k_2)<->(k2) 2 A3", "A3 ->(k3) "])
        >>> net.print_complex_matrix()
             A1  A2 + A3  A2  2A3  A3
        A1 |  1        0   0    0   0  0 |
        A2 |  0        1   1    0   0  0 |
        A3 |  0        1   0    2   1  0 |
        """
        print_matrix(self.complex_matrix, self.species, list(map(str, self.complexes)))


    def print_influence_matrix(self, var = "g_"):
        """Print the influence matrix, with labels for species and reactions."""
        print_matrix(self.influence_matrix(var), self.species, self.reactionids)


    def inspect(self, print_reactions = False, print_matrices = False, invariants = False):
        """Print information on the chemical reaction network.
        If print_reactions = True, print the network reactions.
        If print_matrices = True, print the stoichiometric, incidence
        and kinetic matrix.
        If invariants = True, print the place and transition invariants (pycddlib required)."""

        print("{} species: {}".format(self.n_species, ", ".join(self.species)))
        print("{} complex{}: {}".format(self.n_complexes,
                                        "es" if self.n_complexes > 1 else "",
                                        ", ".join(map(str, self.complexes))))
        if print_reactions:
            print("{} reaction{}:".format(self.n_reactions,
                                          "s" if self.n_reactions > 1 else ""))
            for r in self.reactions: print(r)
        else:
            print("{} reaction{}.".format(self.n_reactions,
                                          "s" if self.n_reactions > 1 else ""))

        source = self.source_species
        if len(source) > 0:
            print("{} source species: {}".format(len(source), ", ".join(source)))
        else: print("No source species.")

        source = self.source_complexes
        if len(source) > 0:
            print("{} source complex{}: {}".format(len(source),
                                                   "es" if len(source) > 1 else "",
                                                   ", ".join(map(str, source))))
        else: print("No source complexes.")

        sink = self.sink_species
        if len(sink) > 0:
            print("{} sink species: {}".format(len(sink), ", ".join(sink)))
        else: print("No sink species.")

        sink = self.sink_complexes
        if len(sink) > 0:
            print("{} sink complex{}: {}".format(len(sink),
                                                 "es" if len(sink) > 1 else "",
                                                 ", ".join(map(str, sink))))
        else: print("No sink complexes.")

        simple = self.simple_complexes
        if len(simple) > 0:
            print("{} simple complex{}: {}".format(len(simple),
                                                   "es" if len(simple) > 1 else "",
                                                   ", ".join(map(str, simple))))
        else: print("No simple complexes.")

        stoich1 = self.stoich_1_species
        if len(stoich1) > 0:
            print("{} stoichiometry 1 species: {}".format(len(stoich1), ", ".join(stoich1)))
        else: print("No stoichiometry 1 species.")

        constant = self.constant_species
        if len(constant) > 0:
            print("{} constant species: {}".format(len(constant), ", ".join(constant)))
        else: print("No constant species.")

        intermediate = self.intermediate_complexes
        if len(intermediate) > 0:
            print("{} intermediate complex{}: {}".format(len(intermediate),
                                                         "es" if len(intermediate) > 1 else "",
                                                         ", ".join(map(str, intermediate))))
        else: print("No intermediate complexes.")

        intermediate_simple = self.intermediate_simple_complexes
        if len(intermediate_simple) > 0:
            print("{} intermediate simple complex{}: {}".format(len(intermediate_simple),
                                                                "es" if len(intermediate_simple) > 1 else "",
                                                                ", ".join(map(str, intermediate_simple))))
        else: print("No intermediate simple complexes.")

        intermediate_stoich1 = self.intermediate_stoich_1_species
        if len(intermediate_stoich1) > 0:
            print("{} intermediate stoichiometry 1 species: {}".format(len(intermediate_stoich1),
                                                                       ", ".join(intermediate_stoich1)))
        else: print("No intermediate stoichiometry 1 species.")

        if len(self.removed_species) > 0:
            print("Removed species:")
            for s, f in self.removed_species: print("{} = {}".format(s,f))

        print("Stoichiometric space dimension: {}".format(self.stoich_space_dim))
        conservations = self.cons_laws
        if len(conservations) == 1:
            print("One conservation law: {}".format(conservations[0]))
        else:
            if len(conservations) > 1: print("Conservation laws: {}".format(", ".join(map(str, conservations))))
            else: print("No conservations.")

        if invariants:
            pinvariants = self.format_p_invariants()
            if len(pinvariants) == 1:
                print("One P-invariant: {}".format(pinvariants[0]))
            else:
                if len(pinvariants) > 1: print("P-invariants: {}".format(", ".join(map(str, pinvariants))))
                else: print("No P-invariants.")
            tinvariants = self.format_t_invariants()
            if len(tinvariants) == 1:
                print("One T-invariant: {}".format(tinvariants[0]))
            else:
                if len(tinvariants) > 1: print("T-invariants: {}".format(", ".join(map(str, tinvariants))))
                else: print("No T-invariants.")

        print("Network deficiency: {}".format(self.deficiency))
        print("Is mass action: {}".format(self.is_ma))

        if print_matrices:
            print("Stoichiometric matrix:")
            self.print_stoich_matrix()
            print("Incidence matrix:")
            self.print_incidence_matrix()
            print("Kinetic matrix:")
            self.print_kinetic_matrix()
        print("Equations:")
        for e in self.format_equations(): print(e)


    ### Output to file ###

    def save_sbml(self, filepath):
        """Save the network as an SBML."""
        self.logger.info("Going to save model to {}".format(filepath))
        if not self.model:
            self._model, self._document, _ = model_from_reacts(self.reactions)
        success = writeSBMLToFile(self.document, filepath)
        if not success:
            raise SystemExit("Error trying to write SBML model to {}.".format(filepath))


    def save_reaction_file(self, filepath, rate = False, precision = 3):
        """Save the reactions to filepath."""
        with open(filepath, 'w') as f:
            for r in self.reactions:
                f.write(r.format(rate, precision) + '\n')


    def save_to_file(self, filepath, reactions = None, overwrite = 'a', rs = False, log = None):
        """Write reactions (with optional removed species) to filepath,
        with an optional line of log."""
        with open(filepath, overwrite) as f:
            if log:
                if reactions:
                    df = from_reacts(reactions).format_deficiency()
                else:
                    df = self.format_deficiency()
                log = log + ", " + df
                f.write(log + '\n')

            if not reactions:
                reactions = self.reactions
            for r in reactions:
                f.write(str(r) + '\n')

            if rs and len(self.removed_species) > 0:
                f.write("Removed species" + '\n')
                for s, e in self.removed_species:
                    f.write(s + " = " + str(e) + '\n')
            f.write('\n')


def from_sbml(xmlfile):
    """Create a CRN from an SBML file."""
    return CRN(model_from_sbml(xmlfile))


def from_react_file(filename, rate = False):
    """Create a CRN from a reaction file."""
    return from_reacts(parse_reaction_file(filename, rate))


def from_reacts(reacts):
    """Create a CRN from a list of reactions."""
    return CRN.from_reacts(reacts)


def from_react_strings(reacts, rate = False):
    """Create a CRN from a list of reaction strings."""
    return from_reacts(parse_reactions(reacts, rate))


def simulate_crn(rates, initials, molecular_weights, end_time=4, crn=None, incr=0.001, v=1, return_mass_fraction=True):
    """Simulate the deterministic dynamics."""
    # checking the conservation of mass in all reactions
    assert_cons_law(crn, molecular_weights)
    times = np.arange(0, end_time, incr)
    par = dict(zip(crn.kinetic_params, rates))
    # inserting rate constants in derivatives
    eqs = [e.subs(par.items()) for e in crn.equations()]
    # turning sympy equations into lambda functions
    lam_funcs = list(map(lambda eq: sp.lambdify(crn.species, eq), eqs))
    # integrate and multiply by v to get moles
    moles = integrate.odeint(lambda x, t: list(map(lambda func: func(*x), lam_funcs)), initials, times) * v
    if return_mass_fraction:
        # convert to mass
        masses = moles * np.array(molecular_weights)[np.newaxis,...]
        mass_fractions = masses / np.sum(masses, axis=1, keepdims=True)
        return times, mass_fractions
    else:
        mol_fractions = moles / np.sum(moles, axis=1, keepdims=True)
        return times, mol_fractions


def assert_cons_law(crn, molecular_weights):
    net_masses = np.array(molecular_weights)[np.newaxis, ...] @ crn.stoich_matrix
    for i, rxn_net_mass in enumerate(net_masses):
        assert rxn_net_mass == 0, f'Mass is not conserved in rxn {i+1}! Check your input crn and species molecular weights.'
