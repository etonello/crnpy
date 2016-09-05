#!/usr/bin/env python

"""CRN class."""

from collections import defaultdict
import itertools
from libsbml import writeSBMLToFile, formulaToL3String, SBMLDocument
import logging
import numpy as np
from scipy.sparse.csgraph import connected_components
import sympy as sp
import warnings

from .conslaw import ConsLaw
from .createmodel import model_from_reacts, model_from_sbml, replace_reacts
from .matrixfunctions import negative, sdiag, print_matrix, _pos_dependent, _pos_generators, dependent
from .crncomplex import Complex, sympify
from .parsereaction import parse_reaction_file, parse_complex, parse_reactions
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
    >>> net = from_react_strings(["a -> b", "2b -> "])
    """

    def __init__(self, model = None):
        self.logger = logging.getLogger("crnpy.crn")
        self.logger.info("Creating an instance of crn.")

        if not model:
            document = SBMLDocument(3, 1)
            model = document.createModel(), document

        self._model, self._document = model

        # species
        self._species_from_sbml()
        self._removed_species = []

        # reactions
        self._get_reactions()
        self._populate()


    @classmethod
    def from_reacts(cls, reacts):
        crn = cls()
        crn._model, crn._document, crn._species = model_from_reacts(reacts)
        crn._n_species = len(crn._species)
        crn._removed_species = []
        crn._reactions = reacts
        crn._populate()
        return crn


    @property
    def model(self):
        return self._model


    @property
    def document(self):
        return self._document


    @property
    def species(self):
        return tuple(self._species)


    @property
    def removed_species(self):
        """Couples (species, expression) for species that have been eliminated."""
        return tuple(self._removed_species)


    @property
    def complexes(self):
        return tuple(self._complexes)


    @property
    def reactions(self):
        return tuple(self._reactions)
    #@reactions.setter
    #def reactions(self, value):
        #self._reactions = value
        #self._populate()


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
        return self._rates


    @property
    def kinetic_params(self):
        return tuple(self._kinetic_params)


    def _species_from_sbml(self):
        """Extract species from SBML model."""
        # species
        self._species = [self.model.getSpecies(s).getName()
                        if self.model.getSpecies(s).getName()
                        else self.model.getSpecies(s).getId() for s in range(self.model.getNumSpecies())]
        self._species_ids = dict(zip([self.model.getSpecies(s).getId()
                                      if self.model.getSpecies(s).getId()
                                      else self.model.getSpecies(s).getName() for s in range(self.model.getNumSpecies())], self._species))
        self._species = sorted(self.species)
        self._n_species = len(self.species)


    def _get_reactions(self):
        """Extract reactions from SBML model."""
        nr = self.model.getNumReactions()
        reactions = []
        self.__reactionids = []

        for r in range(nr):
            reaction = self.model.getReaction(r)
            reactionid = reaction.getName() if reaction.getName() else reaction.getId()
            self.__reactionids.append(reaction.getId() if reaction.getId() else reaction.getName())

            reactant = Complex(dict((self._species_ids[c.getSpecies()], 1 if np.isnan(c.getStoichiometry()) else int(c.getStoichiometry())) \
                               for c in reaction.getListOfReactants()))
            product = Complex(dict((self._species_ids[c.getSpecies()], 1 if np.isnan(c.getStoichiometry()) else int(c.getStoichiometry())) \
                              for c in reaction.getListOfProducts()))

            # remove species with stoichiometric coefficient equal to 0
            reactant = Complex(dict((s, sc) for s, sc in reactant.items() if sc != 0))
            product = Complex(dict((s, sc) for s, sc in product.items() if sc != 0))

            # in first approximation, we are ignoring function definitions
            if reaction.getKineticLaw().getMath():
                kineticlaw = sympify(formulaToL3String(reaction.getKineticLaw().getMath()))
                # replace species id with species name
                for sid in self._species_ids:
                    if sid != self._species_ids[sid]:
                        kineticlaw = kineticlaw.subs(sympify(sid), sympify(self._species_ids[sid]))
                # replace parameter id with parameter name
                params = list(reaction.getKineticLaw().getListOfParameters()) + list(self.model.getListOfParameters())
                for param in params:
                    if param.getId() and param.getName():
                        if param.getId() != param.getName():
                            kineticlaw = kineticlaw.subs(sympify(param.getId()), sympify(param.getName()))

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
                rate = reaction.reactant.ma() * sympify(reaction.getKineticLaw().getListOfParameters().get(0).getId())

            reactions.append(Reaction(reactionid, reactant, product, rate))
            if reaction.getReversible():
                reactions.append(Reaction(reactionid + "_rev", product, reactant, raterev))
                self.__reactionids.append(reaction.getId() + "_rev" if reaction.getId() else reaction.getName() + "_rev")
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


    def replace_reacts(self, new_reactions):
        """Replace reactions with newreactions.
        Update the sbml model and all properties."""
        self._model, self._document, self._species = \
            replace_reacts(self.model, self.document, new_reactions)
        self._reactions = new_reactions
        self._populate()


    def set_params(self, params_dict):
        """Replace the parameters used in the reaction rates
        with the values in dictionary params_dict.

        params_dict is a dictionary with keys the parameters to replace, and
        values the sympy expressions or numeric values that replace them.

        In the following example we set all the kinetic parameters to 0.001:

        :Example:

        >>> from crnpy.crn import CRN, from_react_strings
        >>> net = from_react_strings(["A1 ->(k1) A2 + A3", "A2 ->(k2) 2 A3"])
        >>> net.set_params(dict((k, 1) for k in net.kinetic_params))
        >>> net.reactions
        (r0: A1 ->(1.000e-3) A2 + A3, r1: A2 ->(1.000e-3) 2A3)

        """
        self.replace_reacts([Reaction(r.reactionid,
                                      r.reactant,
                                      r.product,
                                      r.rate.subs(params_dict)) for r in self.reactions])


    ### Matrices ####

    @property
    def complex_matrix(self):
        """Complex matrix."""
        # ns x nc matrix
        return self._complex_matrix


    @property
    def incidence_matrix(self):
        """Incidence matrix."""
        # nc x nr matrix
        return self._incidence_matrix


    @property
    def stoich_matrix(self):
        """Stoichiometric matrix."""
        # s x r matrix
        return sp.SparseMatrix(self.complex_matrix.multiply(self.incidence_matrix))


    @property
    def kinetic_matrix(self):
        """Kinetic matrix."""
        # -Laplacian
        return -self.incidence_matrix.multiply(sdiag(self.kinetic_params).multiply(negative(self.incidence_matrix).T))


    @property
    def laplacian(self):
        """Generalised Laplacian of the graph of complexes."""
        # -kinetic matrix
        return self.incidence_matrix.multiply(sdiag(self.kinetic_params).multiply(negative(self.incidence_matrix).T))


    def influence_matrix(self, var = None, state = None, check = False, interval = None, params = None):
        """Return the n_s x n_r influence matrix.
        The element at position ij is a variable with a plus in front
        if the rate v_j(x) of reaction j increases with x_i,
        a minus if it decreases, and 0 if it is constant in x_i.
        Every reaction rate v_j is assumed to be strictly monotone in each variable x_i.
        Optional check for monotonicity is available with check = True only for
        functions of one variable.

        Keyword arguments:

        * var -- variable name to use. Default is g\_.
        * state -- dictionary of the point where the derivatives are evaluated. Default is (1, 1, ..., 1).
        * check -- check for monotonicity. Only available for unary functions.
        * interval -- interval where monotonicity is checked - default is [0, oo).
        * params -- values for the kinetic parameters. All set to 1 by default.

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
                params = dict((sympify(f), params[f]) for f in params)
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

                deriv = sp.ratsimp(sp.diff(v, sympify(s)))

                # if no state is provided, use (1, 1, ..., 1)
                # and set kinetic parameters to 1
                if not state:
                    state = {}
                state = dict((sympify(f), state[f]) for f in state)
                for f in deriv.free_symbols:
                    if f not in state:
                        state[f] = 1

                # evaluate derivative in state
                deriv = deriv.subs(state)

                if deriv != 0:
                    if deriv > 0:
                        im[i, j] = sympify(var + str(i + 1) + "_" + str(j + 1))
                    else:
                        im[i, j] = sympify("-" + var + str(i + 1) + "_" + str(j + 1))
        return im


    ### Other attributes ###

    @property
    def deficiency(self):
        """Deficiency of the chemical reaction network,
        calculated as number of complexes, minus number of linkage classes,
        minus the rank of the stoichiometric matrix."""
        return sp.Matrix(self.incidence_matrix).rank() - sp.Matrix(self.stoich_matrix).rank()


    @property
    def n_linkage_classes(self):
        """Number of linkage classes of the graph of complexes."""
        return self.n_complexes - sp.Matrix(self.incidence_matrix).rank()


    def format_deficiency(self):
        """Return a string 'deficiency delta = n_c - n_lc - stoich_matrix_rank"""
        return "deficiency {} = {} - {} - {}".format(self.deficiency, \
                                                     self.n_complexes, \
                                                     self.n_linkage_classes, \
                                                     self.stoich_matrix.rank())


    @property
    def is_ma(self):
        """Return True if the network has mass-action kinetics."""
        return all(sympify(s) not in k.cancel().atoms() for s in self.species for k in self.kinetic_params)


    @property
    def stoich_space_dim(self):
        """Return the dimension of the stoichiometric subspace,
        i.e. the rank of the stoichiometric matrix."""
        return sp.Matrix(self.stoich_matrix).rank()


    def _cons_laws(self):
        """Return a base of the nullspace in row echelon form."""
        return [sp.Matrix(row) for row in sp.Matrix([list(c) for c in self.stoich_matrix.T.nullspace()]).rref()[0].tolist()]


    @property
    def cons_laws(self):
        """Tuple of generators of the network conservation laws."""
        return tuple((v.T * sp.Matrix(list(map(sympify, self.species))))[0, 0] for v in self._cons_laws())


    @property
    def p_invariants(self):
        """Return a list of generators of the Petri Net P-invariants.
        Requires pycddlib."""
        return _pos_generators(self.stoich_matrix.T)


    def format_p_invariants(self):
        """Return a list of generators of the Petri Net P-invariants,
        as combinations of the species.
        Requires pycddlib."""
        pinv = self.p_invariants
        if len(pinv) > 0:
            return (pinv * sp.Matrix(list(map(sympify, self.species)))).tolist()
        else:
            return []


    @property
    def t_invariants(self):
        """Return a matrix with rows the generators of the Petri Net T-invariants.
        Requires pycddlib."""
        return _pos_generators(self.stoich_matrix)


    def format_t_invariants(self):
        """Return a matrix with rows the generators of the Petri Net T-invariants,
        as combinations of the reaction ids.
        Requires pycddlib."""
        tinv = self.t_invariants
        if len(tinv) > 0:
            return tinv * sp.Matrix(list(map(sympify, self.reactionids)))
        else:
            return sp.Matrix()


    @property
    def elem_modes(self):
        """Return the list of elementary modes of the network,
        as vectors of length nr."""
        return self.t_invariants.tolist()


    def format_elem_modes(self):
        """Return the list elementary modes of the network,
        as a string with the reaction ids."""
        return self.format_t_invariants().tolist()


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


    def tree_constant(self, index):
        """Return the constant of Prop. 3 of
        ''Toric Dynamical Systems'' by Shul et al.,
        associated to index.
        The term ''tree constant'' is introduced in
        ''Translated Chemical Reaction Networks'' by M. D. Johnston."""
        # find the linkage classes
        _, lcs = self.weak_conn_components()
        l = lcs[index]
        inds = [i for i in range(self.n_complexes) if lcs[i] == l and i != index]
        return (-1)**(self.n_complexes - 1)*self.kinetic_matrix.extract(inds, inds).det()


    def tree_constants(self):
        """Return the constants of Prop. 3 of
        ''Toric Dynamical Systems'' by Shul et al.
        The term ''tree constant'' is introduced in
        ''Translated Chemical Reaction Networks'' by M. D. Johnston."""
        kinetic_matrix = self.kinetic_matrix
        # find the linkage classes
        _, lcs = self.weak_conn_components()
        tree_consts = []
        for index in range(self.n_complexes):
            l = lcs[index]
            inds = [i for i in range(self.n_complexes) if lcs[i] == l and i != index]
            tree_consts.append((-1)**(self.n_complexes - 1)*kinetic_matrix.extract(inds, inds).det())
        return tree_consts


    ### System of ODEs ###

    def equations(self):
        """Return the derivatives of the species concentrations."""
        equations = self.stoich_matrix * self.rates
        return equations


    def print_equations(self):
        """Print the system of ODEs."""
        odes = self.equations()
        for s in range(self.n_species):
            print ("d{}/dt = {}".format(self.species[s], odes[s, 0]))


    def groebner(self):
        """Return a groebner basis for the steady state ideal."""
        return sp.groebner([sp.ratsimp(e).as_numer_denom()[0] for e in self.equations()],
                           *[sympify(s) for s in self.species])


    def _constant_species(self):
        """Constant species indices."""
        # This is more general than
        # [i for i in range(self.n_species) if all([j == 0 for j in self.stoich_matrix[i,:]])]
        eqs = self.equations()
        return [i for i in range(self.n_species) if eqs[i] == 0]


    @property
    def constant_species(self):
        """Constant species."""
        return [self.species[s] for s in self._constant_species()]


    def derivative(self, cplx):
        """Return the derivative expression of a complex."""
        c = parse_complex(cplx)
        der = sympify(0)
        # sum the derivative of each species, multiplied by the stoichiometry
        for s in c:
            if str(s) not in self.species:
                raise ValueError("Invalid complex: {} not a valid species.".format(s))
            else:
                der = der + c[s] * (self.stoich_matrix[self.species.index(str(s)),:] * self.rates)[0,0]
        return der.cancel()


    def is_constant(self, cplx):
        """Return True if the derivative is zero."""
        return self.derivative(cplx) == 0


    def has_linear_equation(self, species):
        """Check if the equation ds/dt = 0 is linear in s."""
        expr = sp.ratsimp((self.stoich_matrix[self.species.index(species), :] * self.rates)[0]).as_numer_denom()[0]
        return sp.degree(expr, sympify(species)) == 1


    def is_dyn_eq(self, net):
        """Check if two networks have the same dynamics.
        They must have the same set of species."""
        if set(self.species) != set(net.species):
            return False
        eqs = self.equations()
        othereqs = net.equations()
        return all((eqs[i] - othereqs[net.species.index(self.species[i])]).cancel() == 0 for i in range(self.n_species))


    def is_emul(self, net, morphism = None):
        """Check if the network is an emulation of the second
        under the provided morphism, i.e. the networks have the
        same equations modulo a renaming of the variables.
        See also Cardelli, L. (2014)., Morphisms of reaction networks
        that couple structure to function. BMC systems biology, 8(1), 84.

        morphim = dictionary species -> species. Default = identity.
        """
        if morphism == None:
            morphism = {}
            for s in net.species: morphism[s] = s

        # check if morphism maps species of net to species
        if not (set(morphism.keys()) == set(net.species) and
                set(morphism.values()).issubset(set(self.species))):
            raise ValueError("Morphism does not map all species of net to species of crn.")

        ode = self.equations()
        equations = {}
        for i in range(self.n_species):
            equations[self.species[i]] = ode[i]

        net_ode = net.equations()
        for s in morphism:
            net_ode = net_ode.subs(sympify(s), sympify(morphism[s]))
        newequations = {}
        for i in range(net.n_species):
            if morphism[net.species[i]] in newequations:
                if net_ode[i] != newequations[morphism[net.species[i]]]:
                    return False
            else: newequations[morphism[net.species[i]]] = net_ode[i]
        return all([(newequations[s] - equations[s]).cancel() == 0 for s in newequations])


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
        """Check whether the species in ss are intermediates,
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

        Keyword arguments:

        * rapid_eq -- list of couples (species to remove, complex to replace species with).
        * qss -- list of species that to remove via qss.
        * cons_law -- couple (string, ConsLaw), where the first element is the species
                      that will be written as a function of the other species in the conservation.
                      If the reduction method is not specified for a species in the conservation,
                      the quasi-steady state approximation is used.
        * minimal -- find network of minimal structure when applying qss.
        * remove_const -- remove any species left constant after the reduction.
        * merge_reacts -- merge reactions with the same reactant and product.
        * adjust -- change the rates so that they reflect the reactant species.
        * network_file -- save the reduction steps to file with given path.

        Update remove_species.
        """
        # Algorithms from Tonello et al. (2016),
        # On the elimination of intermediate species in chemical reaction networks.

        if cons_law:
            add_species = [str(c) for c in cons_law[1].species.keys() if c != sympify(cons_law[0])]
        else: add_species = []

        if rapid_eq != None: rapid_eq_species = [couple[0] for couple in rapid_eq]
        else: rapid_eq_species, rapid_eq = [], []
        if qss != None: qss_species = qss + [a for a in add_species if a not in rapid_eq_species and a not in qss]
        else: qss_species = [a for a in add_species if a not in rapid_eq_species]

        if network_file: self.save_to_file(network_file, overwrite = 'w', log = "Original network")

        # rapid-equilibrium
        for couple in rapid_eq:
            self._rapid_eq(couple, debug)
            if network_file: self.save_to_file(network_file, log = "Rapid eq. on {}, {}".format(couple[0], couple[1]))
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
        if adjust: self._fix_ma(debug)

        if network_file: self.save_to_file(network_file, rs = True, log = "Final network")


    def remove_by_cons(self, species, cons_law, debug = False):
        """Remove a species using a conservation law.
        First replace removed_species in the conservation law with their expression.
        Then use the conservation expression to write the species
        concentration as function of the remaining species.
        """
        # Algorithm from Tonello et al. (2016),
        # On the elimination of intermediate species in chemical reaction networks.
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
        #conservation = (conservation / sympify(species)).cancel()
        #exp = cons_law.constant / conservation
        exp = sp.solve(conservation - cons_law.constant, sympify(species))[0]

        # remove species
        self.remove_constant(species, expr = exp)

        if debug: print("Remove by Conservation: added to removed_species {}".format(self.removed_species))


    def qss(self, intermediate = None, cons_law = None, minimal = False, remove_const = False, \
            merge_reacts = False, adjust = False, debug = False, network_file = None):
        """Remove an intermediate species via quasi-steady state approximation.

        Keyword arguments:

        * intermediate -- species to eliminate.
        * cons_law -- couple (string, ConsLaw), where the first element is the species
                      that will be written as a function of the other species in the conservation.
                      All other species in the conservation are eliminated via quasi-steady state approximation.
        * minimal -- find network of minimal structure when applying qss.
        * remove_const -- remove any species left constant after the reduction.
        * merge_reacts -- merge reactions with the same reactant and product.
        * adjust -- change the rates so that they reflect the reactant species.
        * network_file -- save the reduction steps to file with given path.

        Update remove_species.
        """
        # Algorithm from Tonello et al. (2016),
        # On the elimination of intermediate species in chemical reaction networks.
        # Algorithm for minimal structure based on Madelaine et al. (2026),
        # Normalizing Chemical Reaction Networks by Confluent Structural Simplification. (CMSB 2016).
        return self.remove(qss = ([intermediate] if intermediate else None), cons_law = cons_law,
                           adjust = adjust, minimal = minimal, debug = debug, \
                           network_file = network_file, remove_const = remove_const, merge_reacts = merge_reacts)


    def _qss(self, intermediates, minimal = False, error_if_missing = True, network_file = None, debug = False):
        """Eliminate the intermediates via quasi-steady state approximation."""
        # Algorithm from Tonello et al. (2016),
        # On the elimination of intermediate species in chemical reaction networks.
        # Algorithm for minimal structure based on Madelaine et al. (2026),
        # Normalizing Chemical Reaction Networks by Confluent Structural Simplification. (CMSB 2016).

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

            y = sympify(intermediate)
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
                                                        sp.nsimplify(sympify(coeffs[i]), rational = True)).cancel()
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

        self.replace_reacts(reactions)


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

        y = sympify(intermediate)
        expr = y
        # if there are no reactant of product reactions, species is not an intermediate
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
                                                     sympify("k_" + newid) * newreactant.ma() if no_rates
                                                     else (sp.gcd(n, m) * r1.rate * r2.rate / denom).cancel()))

            for r in range(len(newreactions)):
                newreactions[r]._rate = (newreactions[r].rate).subs(y, expr).cancel()

        if debug:
            print("New reactions")
            for r in newreactions: print(r)
            print

        self.replace_reacts(newreactions)
        self._removed_species.append((intermediate, expr))


    def _fix_ma(self):
        """Adjust the reactant and products so that they match the
        monomials in the numerators of the reaction rates."""
        reactions = list(itertools.chain(*map(lambda r: _split_reaction_monom(r, self.species), self.reactions)))
        for r in reactions: r._fix_ma(self.species)
        self.replace_reacts(reactions)


    def _remove_react_prod(self):
        """For each reaction, remove the species that
        are in common between reactant and product."""
        for r in self.reactions: r.remove_react_prod()
        self.replace_reacts(self.reactions)


    def _same_denom(self):
        self.replace_reacts(_same_denom(self.reactions))


    def _fix_denom(self):
        """Remove reactant and product species if their
        concentrations divides the denominator of the rate."""
        reactions = list(itertools.chain(*map(lambda r: _split_reaction_monom(r, self.species), self.reactions)))
        for r in reactions: r._fix_denom(self.species)
        self.replace_reacts(reactions)


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
        keep_pool = sympify(pool_name) / (1 + exp / sympify(str(keep)))
        keep_factor = (1 / (1 + exp / sympify(str(keep)))).cancel()
        remove_pool = keep_pool * exp / sympify(str(keep))
        remove_factor = (keep_factor * exp / sympify(str(keep))).cancel()

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
                rate = (keep_factor * rate * pool_complex.ma() / sympify(str(keep))).cancel()
            if str(reaction.reactant) == remove:
                reactant = pool_complex
                rate = (remove_factor * rate * pool_complex.ma() / sympify(remove)).cancel()
            if reactant != product:
                reactions.append(Reaction(reaction.reactionid, \
                                          reactant, \
                                          product, \
                                          rate.cancel()))

        self.replace_reacts(reactions)
        self._removed_species.append((remove, exp))

        if debug:
            for r in self.reactions: print(r)

        if cons_law:
            self.remove_by_cons(cons_law[0], cons_law[1], debug)


    def rapid_eq(self, recouple, cons_law = None, debug = False, network_file = None):
        """Apply the rapid equilibrium approximation to recouple,
        replacing the species in recouple[0] with the complex in recouple[1].

        Keyword arguments:

        * recouple -- couple (species to eliminate, complex to replace species)
        * cons_law -- couple (string, ConsLaw), where the first element is the species
                      that will be written as a function of the other species in the conservation.
                      All other species in the conservation are eliminated via quasi-steady state approximation.
        * network_file -- save the reduction steps to file with given path.

        Update remove_species.
        """
        # Algorithm from Tonello et al. (2016),
        # On the elimination of intermediate species in chemical reaction networks.
        return self.remove(rapid_eq = [recouple], cons_law = cons_law, debug = debug, network_file = network_file)


    def _rapid_eq(self, couple, debug):
        """The reactions between the complexes in couple
        are assumed at rapid equilibrium.
        The first complex is removed and the second
        replaces it in the network."""
        # Algorithm from Tonello et al. (2016),
        # On the elimination of intermediate species in chemical reaction networks.

        remove, keep = couple
        removecomplex, keepcomplex = map(parse_complex, couple)

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
            newreactions[r]._rate = (newreactions[r].rate).subs(sympify(remove), exp)

        self.replace_reacts(newreactions)

        self._removed_species.append((remove, exp))
        if debug: print("Rapid Eq: added to removed_species {}".format(self.removed_species))


    def merge_network(self, network):
        """Add reactions of network."""
        if len(network.n_reactions) > 0:
            self.replace_reacts(self.reactions + network.reactions)


    def detach_network(self, inter_complexes):
        """Take a list of complexes, and remove the weakly connected components
        for these complexes from the network, returning them as a separate network."""
        connComponents = self.weak_conn_components()[1]
        involvedComponents = set([connComponents[i] for i in inter_complexes])
        involvedComplexes = [i for i in range(self.n_complexes) if connComponents[i] in involvedComponents]

        involvedReactions = [r for r in range(self.n_reactions) if list(self.incidence_matrix.col(r)).index(-1) in involvedComplexes \
                                                                      or list(self.incidence_matrix.col(r)).index(1) in involvedComplexes]

        if len(involvedComplexes) != self.n_complexes:
            allreactions = self.reactions

            network = CRN.from_reacts([allreactions[r] for r in involvedReactions])

            reactions = [allreactions[r] for r in range(self.n_reactions) if not r in involvedReactions]
            self.replace_reacts(reactions)
            return network
        return None


    def remove_constant(self, const_species, expr = None, debug = False):
        """Remove the species from the network by setting it to constant.
        Give a warning if the species is not constant in the original network.

        * expr -- optional expression that replaces const_species in rates.

        """
        # Algorithm from Tonello et al. (2016),
        # On the elimination of intermediate species in chemical reaction networks.

        # check that const_species is a valid species
        # either const_species is in set of species, and has null derivative
        # or is used in kinetics
        if const_species not in self._species + [str(v) for v in list(set(list(itertools.chain(*[k.atoms() for k in self.kinetic_params]))))]:
            warnings.warn("Species {} not valid.".format(const_species))
        else:
            if const_species in self._species and self.derivative(const_species) != 0:
                warnings.warn("Species {} not constant.".format(const_species))

        self.logger.info("Removing constants: {}".format(', '.join(const_species)))

        # remove species from reactants and products
        # and replace species with expression in rates
        for reaction in self.reactions:
            if const_species in reaction.reactant:
                del reaction.reactant[const_species]
                reaction._rate = reaction.rate
            if const_species in reaction.product:
                del reaction.product[const_species]
            if expr != None: reaction._rate = reaction.rate.subs(sympify(const_species), expr).cancel()

        reactions = [r for r in self.reactions if r.reactant != r.product]
        self.replace_reacts(reactions)
        self._removed_species = self._removed_species + [(const_species, expr if expr else sympify(const_species))]


    def remove_all_constants(self, debug = False):
        """Remove all constant species."""
        for x in self.constant_species:
            self.remove_constant(x, debug = debug)


    def merge_reactions(self):
        """Merge the reactions with same reactant and same product."""
        self.replace_reacts(merge_reactions(self.reactions))


    ### Graphs ###

    def complex_graph_adj(self):
        """Return the nc x nc adjacency matrix of the graph of complexes."""
        adjacency = - self.incidence_matrix.multiply(negative(self.incidence_matrix).T)
        for i in range(self.n_complexes): adjacency[i, i] = 0
        adjacency = np.ma.masked_values(adjacency, 0)
        return adjacency


    def dsr_graph_adj(self, keep = None):
        """Return the (ns + nr) x (ns + nr) adjacency matrix
        of the directed species-reaction graph.
        Optionally keep only variables in keep."""
        A = self.stoich_matrix
        im = self.influence_matrix()
        if keep is not None:
            im = im.applyfunc(lambda x: x if x in keep else 0)
        im = im.applyfunc(lambda x: -1 if str(x)[0] == "-" else (1 if x != 0 else 0))
        adjacency = sp.zeros(im.rows, im.rows).row_join(im).col_join(A.T.row_join(sp.zeros(A.cols, A.cols)))
        adjacency = np.array(adjacency.tolist()).astype(np.float64)
        return np.array(adjacency)


    ### Connectivity properties ###

    def strong_conn_components(self):
        """Return the number of strongly connected components of the graph of complexes,
        and an array of labels for the strongly connected components."""
        adjacency = self.complex_graph_adj()
        return connected_components(adjacency, connection = 'strong')


    def strong_terminal_conn_components(self):
        """Return the list labels for the terminal strongly connected components,
        of the graph of complexes, and the array of labels for the strongly connected components."""
        n, cs = self.strong_conn_components()
        adjacency = self.complex_graph_adj()
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


    def terminal_complexes(self):
        """Return the list of complexes in the terminal strongly connected components."""
        return [self.complexes[c] for c in self._terminal_complexes()]


    def _non_terminal_complexes(self):
        """Return the indices of the complexes in the
        nonterminal strongly connected components."""
        stcc, cs = self.strong_terminal_conn_components()
        return [c for c in range(self.n_complexes) if cs[c] not in stcc]


    def non_terminal_complexes(self):
        """Return the list of complexes in the nonterminal strongly connected components."""
        return [self.complexes[c] for c in self._non_terminal_complexes()]


    def weak_conn_components(self):
        """Return the number of weakly connected components of the graph of complexes,
        and an array of labels for the weakly connected components."""
        return connected_components(self.complex_graph_adj())


    def linkage_classes(self):
        """List of complexes grouped by linkage class."""
        nlc, lcs = self.weak_conn_components()
        return [[self.complexes[c] for c in range(self.n_complexes) if lcs[c] == lc] for lc in range(nlc)]


    @property
    def is_weakly_rev(self):
        """Return True if the network is weakly reversible."""
        return self.strong_conn_components()[0] == self.weak_conn_components()[0]


    ### Robustness ###

    def acr_species(self, subnets = False):
        """Return some of the species with absolute concentration robustness.

        * subnets -- if True, look for a decomposition in subnetworks using split_by_ems,
                   and use the decomposition to look for robust ratios and species.

        """
        # Algorithm based on the presentation in the supplementary material of
        # Shinar, G., Feinberg, M., Structural sources of robustness in biochemical reaction networks. Science 2010.
        acr_s = []
        acr_c = self.acr_complexes(as_vectors = True, subnets = subnets)
        if acr_c:
            for i in range(self.n_species):
                if dependent(acr_c, [0 if j != i else 1 for j in range(self.n_species)]):
                    acr_s.append(self.species[i])
        return acr_s


    def acr_complexes(self, as_vectors = False, subnets = False):
        """Return complex ratios that are robust.

        * as_vectors -- if True, return vectors of length = n_species corresponding to the ratios.
        * subnets -- if True, look for a decomposition in subnetworks using split_by_ems,
                     and use the decomposition to look for robust ratios.

        """
        # Algorithm based on the presentation in the supplementary material of
        # Shinar, G., Feinberg, M., Structural sources of robustness in biochemical reaction networks. Science 2010.
        if self.is_ma:
            if subnets:
                nets = self.subnets()
                acr_cs = []
                for net in nets:
                    if as_vectors:
                        def convert(v, ss):
                            # converts vector in subnetwork to vector in crn
                            return [v[ss.index(self.species[i])] if self.species[i] in ss else 0 for i in range(len(self.species))]
                        acr_cs = acr_cs + list(map(lambda v: convert(v, net.species), net.acr_complexes(as_vectors = as_vectors)))
                    else:
                        acr_cs = acr_cs + net.acr_complexes(as_vectors = as_vectors)
                return acr_cs
            else:
                if self.deficiency == 1:
                    ntcs_inds = self._non_terminal_complexes()
                    if as_vectors:
                        return list(set([tuple(self.complexes[i].to_vector(self.species) -
                                               self.complexes[j].to_vector(self.species)) for i in ntcs_inds for j in ntcs_inds if j > i]))
                    return list(set([self.complexes[i].ma()/self.complexes[j].ma() for i in ntcs_inds for j in ntcs_inds if j > i]))
                if self.deficiency == 0 and self.is_weakly_rev:
                    nlc, lc = self.weak_conn_components()
                    ms = []
                    for l in range(nlc):
                        inds = [i for i in range(self.n_complexes) if lc[i] == l]
                        if as_vectors:
                            ms = ms + [tuple(self.complexes[i].to_vector(self.species) -
                                             self.complexes[j].to_vector(self.species)) for i in inds for j in inds if j > i]
                        else:
                            ms = ms + [self.complexes[i].ma()/self.complexes[j].ma() for i in inds for j in inds if j > i]
                    return list(set(ms))
        return []


    ### Subnetworks ###

    def subnets(self):
        """ Try to split the network into subnetworks such that
        the ranks of the subnetworks sum to the rank of the network.
        Use split_by_ems to look for candidate subnetworks."""
        reacts, rem_reacts = self.split_by_ems()
        nets = [from_reacts(r) for r in reacts] + ([from_reacts(rem_reacts)] if rem_reacts else [])
        if sum([net.stoich_matrix.rank() for net in nets]) == self.stoich_matrix.rank():
            return nets
        else:
            return [self]


    def split_by_ems(self, same_react = False):
        """ Split reactions according to the elementary modes.
        If same_react = True, do not split reactions with the same reactant."""
        subnet = {}
        reacts = {}
        tinvs = self.t_invariants

        if tinvs.rows == 0:
            print("No elementary modes.")
            return [self.reactions], None

        for r in range(self.n_reactions):
            if tinvs[0, r] != 0:
                subnet[r] = 0
                reacts[str(self.reactions[r].reactant)] = 0
        for t in range(1, tinvs.rows):
            inds = [r for r in range(self.n_reactions) if tinvs[t, r] != 0]

            values = [subnet[h] for h in inds if h in subnet]
            if same_react:
                values = values + [reacts[str(self.reactions[h].reactant)] for h in inds
                                   if str(self.reactions[h].reactant) in reacts]

            if len(values) > 0: value = min(values)
            else: value = max(subnet.values()) + 1

            for h in inds:
                subnet[h] = value
                reacts[str(self.reactions[h].reactant)] = value

        # case of reactions not taking part in any em
        reacts_not_ems = []
        for r in range(self.n_reactions):
            if r not in subnet:
                warnings.warn("Network not covered by elementary modes.")
                reacts_not_ems.append(self.reactions[r])

        return [[self.reactions[r] for r in [h for h in range(self.n_reactions) if h in subnet and subnet[h] == sb]]
                    for sb in range(max(subnet.values()) + 1)], reacts_not_ems if len(reacts_not_ems) > 0 else None


    ### Print ###

    def print_laplacian(self, numeric = False, precision = 3):
        """Print the laplacian."""
        print_matrix(self.laplacian,
                     list(map(str, self.complexes)),
                     list(map(str, self.complexes)),
                     numeric,
                     precision)

    def print_kinetic_matrix(self, numeric = False, precision = 3):
        """Print the laplacian."""
        print_matrix(self.kinetic_matrix,
                     list(map(str, self.complexes)),
                     list(map(str, self.complexes)),
                     numeric,
                     precision)


    def print_stoich_matrix(self):
        """Print the stoichiometrix matrix."""
        print_matrix(self.stoich_matrix, self.species, self.reactionids)


    def print_incidence_matrix(self):
        """Print the incidence matrix."""
        print_matrix(self.incidence_matrix, list(map(str, self.complexes)), self.reactionids)


    def print_complex_matrix(self):
        """Print the complex matrix."""
        print_matrix(self.complex_matrix, self.species, list(map(str, self.complexes)))


    def print_influence_matrix(self):
        """Print the influence matrix."""
        print_matrix(self.influence_matrix(), self.species, self.reactionids)


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
        self.print_equations()


    ### Output to file ###

    def save_sbml(self, filepath):
        """Save the network as an SBML."""
        self.logger.info("Going to save model to {}".format(filepath))
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
