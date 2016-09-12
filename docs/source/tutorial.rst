crnpy tutorial
==============

Creating a chemical reaction network
------------------------------------

A chemical reaction network (CRN) object can be created from an SBML
file:

.. code:: python

    >>> from crnpy.crn import CRN, from_sbml, from_react_strings, from_react_file
    >>> enz = from_sbml("examples/data/sbml/enzyme.xml")

or by specifying a list of strings describing reactions in a human-readable format:

.. code:: python

    >>> reactions = ['e + a (k_0)<->(k0) ea',
    ...              'e + b (k_1)<->(k1) eb',
    ...              'ea + b (k_2)<->(k2) eab',
    ...              'eb + a (k_3)<->(k3) eab',
    ...              'eab -> e + p']
    >>> bi_uni_random = from_react_strings(reactions)

The list of reactions can be read from a text file, using
*from\_react\_file*:

.. code:: python

    bio26 = from_react_file("examples/data/reactions/biomodels/biomd0000000026")

The folder *examples/data/reactions/* contains examples of networks in human-readable format,
with the subfolder *biomodels* containing examples
from the `BioModels <http://biomodels.caltech.edu/>`_ database [7]_.

A string in this format contains the reactant and product complexes,
separated either by the characters "->"
for unidirectional reactions, or by the characters "<->" for reversible
reactions.
Notice however that reversible reactions will always be converted to two separate
reactions.

The symbol ">" can be followed by an expression in parenthesis, which
will be, by default, interpreted as a generalised kinetic parameter. In
other words, the reaction rate will be set to the expression in
parenthesis multiplied by the concentrations of the reactant species. In
alternative, one can specify the rate between parenthesis, and use
*from\_react\_strings* or *from\_react\_file* with the option *rate =
True*. For example, the network created with the command

.. code:: python

    >>> one_react_net = from_react_strings(["a + b ->(k1) c"])

coincides with the network created by the command

.. code:: python

    >>> one_react_net = from_react_strings(["a + b ->(k1*a*b) c"], rate = True)

The expression in parenthesis must be a string that can be successfully
converted to a *SymPy* expression. Numeric values are of course
accepted. Stoichiometric coefficients in the reactant and product
complexes are instead specified by integers preceding the species, as in
the following example:

.. code:: python

    >>> one_react_net = from_react_strings(["2a + b ->(k*a**2*b/(a+b+c)) 3 c"], rate = True)

The symbol "<" can also be preceded by an expression in parenthesis,
which will be used to set the rate of the reverse reaction.

A reaction id can be specified before the reactant, and must be followed
by a colon. For example, we can specify an id for the last reaction with

.. code:: python

    >>> one_react_net = from_react_strings(["r_in: 2a + b ->(k*a**2*b/(a+b+c)) 3 c"], rate = True)

If an id is not specified, reactions are assigned an id of the form
*r\_n*, with n an integer, starting from 0.
For example, for the bi_uni_random example above we have:

.. code:: python

    >>> for r in bi_uni_random.reactions: print(r)
    ... 
    r0: a + e ->(k0) ea
    r0_rev: ea ->(k_0) a + e
    r1: b + e ->(k1) eb
    r1_rev: eb ->(k_1) b + e
    r2: b + ea ->(k2) eab
    r2_rev: eab ->(k_2) b + ea
    r3: a + eb ->(k3) eab
    r3_rev: eab ->(k_3) a + eb
    r4: eab ->(k_r4) e + p

Notice that for each reversible reaction two separate reactions have been created,
with *_rev* being appendend to the reaction id to form the id of the reverse reaction.

As shown in the last example, kinetic parameters are optional:
since no kinetic parameter was specified for the last reaction,
it was assigned a parameter equal to *k_* followed by the reaction id.

Comments can be added to a reaction file using the symbol "#". Anything
appearing after the hash sign will be ignored.

Exploring chemical reaction networks
------------------------------------

The library implements some elements of chemical reaction network theory
(see for example [1]_, [3]_, [5]_).

Species, complexes and reactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Attributes of a CRN object include the network species, complexes, and
reactions:

.. code:: python

    >>> enz.species, enz.complexes
    (('E', 'ES', 'P', 'S'), (E + S, ES, E + P))
    >>> for r in enz.reactions: print(r)
    ... 
    veq: E + S ->(comp*veq_kon) ES
    veq_rev: ES ->(comp*veq_koff) E + S
    vcat: ES ->(comp*vcat_kcat) E + P
    >>> bio26.species
    ('M', 'MAPKK', 'MKP', 'M_MAPKK', 'M_MKP', 'Mp', 'Mp_MAPKK', 'Mp_MKP', 'Mp_MKP_', 'Mpp', 'Mpp_MKP')
    >>> bio26.complexes
    (M + MAPKK, M_MAPKK, MAPKK + Mp, Mp_MAPKK, MAPKK + Mpp, MKP + Mpp, Mpp_MKP, Mp_MKP, MKP + Mp, Mp_MKP_, M_MKP, M + MKP)
    >>> for r in bio26.reactions: print(r)
    ... 
    binding_MAPK_and_PP_MAPKK: M + MAPKK ->(k1*uVol) M_MAPKK
    binding_MAPK_and_PP_MAPKK_rev: M_MAPKK ->(k_1*uVol) M + MAPKK
    phosphorylation_of_MAPK: M_MAPKK ->(k2*uVol) MAPKK + Mp
    binding_PP_MAPKK_and_P_MAPK: MAPKK + Mp ->(k3*uVol) Mp_MAPKK
    binding_PP_MAPKK_and_P_MAPK_rev: Mp_MAPKK ->(k_3*uVol) MAPKK + Mp
    phosphorylation_of_P_MAPK: Mp_MAPKK ->(k4*uVol) MAPKK + Mpp
    binding_MKP_and_PP_MAPK: MKP + Mpp ->(h1*uVol) Mpp_MKP
    binding_MKP_and_PP_MAPK_rev: Mpp_MKP ->(h_1*uVol) MKP + Mpp
    dephosphorylation_of_PP_MAPK: Mpp_MKP ->(h2*uVol) Mp_MKP
    dissociation_of_MKP_from_P_MAPK: Mp_MKP ->(h3) MKP + Mp
    dissociation_of_MKP_from_P_MAPK_rev: MKP + Mp ->(h_3) Mp_MKP
    binding_MKP_and_P_MAPK: MKP + Mp ->(h4*uVol) Mp_MKP_
    binding_MKP_and_P_MAPK_rev: Mp_MKP_ ->(h_4*uVol) MKP + Mp
    dephosphorylation_of_P_MAPK: Mp_MKP_ ->(h5*uVol) M_MKP
    dissociation_of_MKP_from_MAPK: M_MKP ->(h6*uVol) M + MKP
    dissociation_of_MKP_from_MAPK_rev: M + MKP ->(h_6*uVol) M_MKP

While the species are simple strings, the complexes and reactions are special objects of type
*Complex* and *Reaction* respectively.

An object of type *Complex* is a *Counter*, a python dictionary where the keys are
the species, and the values are the stoichiometric coefficients of the species
in the complex.
Therefore, a complex E + 2S can be defined in crnpy for example as

.. code:: python

    >>> c = Complex({'E': 1, 'S': 2})
    >>> c
    E + 2S

or more briefly with

.. code:: python

    >>> c = Complex(E=1, S=2)

A *Reaction* object can be created by specifying a reaction id, a reactant complex,
a product complex and the reaction rate. The rate must be a SymPy expression,
or a string that can be successfully converted to a SymPy expression:

.. code:: python

    >>> r = Reaction('r_1', Complex(E=1, S=1), Complex(C=1), "k1*E*S")

One can access for example the reactant, product and rate of the reaction:

.. code:: python

    >>> r.reactant
    E + S
    >>> r.product
    C
    >>> r.rate
    E*S*k1

Network matrices
~~~~~~~~~~~~~~~~

Attributes are available to create the main matrices associated to the reaction network.
Available matrices are the
stoichiometric matrix *stoich\_matrix*, the matrix of stoichiometric
coefficients *complex\_matrix* (often called Y in the literature),
the incidence matrix of the complex graph *incidence\_matrix*.
the Laplacian of the graph of complexes *laplacian*, and its negation *kinetic_matrix*.

.. code:: python

    >>> enz.stoich_matrix()
    Matrix([
    [-1,  1,  1],
    [ 1, -1, -1],
    [ 0,  0,  1],
    [-1,  1,  0]])
    >>> enz.complex_matrix
    Matrix([
    [1, 0, 1],
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 0]])
    >>> enz.incidence_matrix
    Matrix([
    [-1,  1,  0],
    [ 1, -1, -1],
    [ 0,  0,  1]])
    >>> enz.laplacian
    Matrix([
    [ _comp*veq_kon,                  -_comp*veq_koff, 0],
    [-_comp*veq_kon, _comp*vcat_kcat + _comp*veq_koff, 0],
    [             0,                 -_comp*vcat_kcat, 0]])


Special methods are available to print some matrices. For example, for
the stoichiometry matrix and the laplacian:

.. code:: python

    >>> enz.print_stoich_matrix()
         veq  veq_rev  vcat
    E  |  -1        1     1 |
    ES |   1       -1    -1 |
    P  |   0        0     1 |
    S  |  -1        1     0 |
    >>> enz.print_laplacian()
                     E + S                                ES  E + P
    E + S |  _comp*veq_kon                   -_comp*veq_koff      0 |
    ES    | -_comp*veq_kon  _comp*vcat_kcat + _comp*veq_koff      0 |
    E + P |              0                  -_comp*vcat_kcat      0 |

Network dynamics
~~~~~~~~~~~~~~~~

The method *odes()* returns the SymPy differential equations describing the evolution of
species concentrations. These can be printed more nicely with *format_equations()*:

.. code:: python

	>>> for eq in enz.odes(): print(eq)
	...
	Eq(Derivative(E(t), t), _comp*vcat_kcat*ES(t) + _comp*veq_koff*ES(t) - _comp*veq_kon*E(t)*S(t))
	Eq(Derivative(ES(t), t), -_comp*vcat_kcat*ES(t) - _comp*veq_koff*ES(t) + _comp*veq_kon*E(t)*S(t))
	Eq(Derivative(P(t), t), _comp*vcat_kcat*ES(t))
	Eq(Derivative(S(t), t), _comp*veq_koff*ES(t) - _comp*veq_kon*E(t)*S(t))
	>>> for e in enz.format_equations(): print(e)
	... 
	dE/dt = -E*S*_comp*veq_kon + ES*_comp*vcat_kcat + ES*_comp*veq_koff
	dES/dt = E*S*_comp*veq_kon - ES*_comp*vcat_kcat - ES*_comp*veq_koff
	dP/dt = ES*_comp*vcat_kcat
	dS/dt = -E*S*_comp*veq_kon + ES*_comp*veq_koff

One can also look at the conservation laws:

.. code:: python

    >>> enz.cons_laws
    (E - P - S, ES + P + S)

or check if two networks are dynamically equivalent:

.. code:: python

    >>> net1 = from_react_strings(['a ->(k) a + 2b'])
    >>> net2 = from_react_strings(['a ->(2*k) a + b'])
    >>> net1.is_dyn_eq(net2)
    True

We can look for a Gröbner basis for the steady state ideal with the method *groebner()*:

.. code:: python

    >>> bio26.groebner()
    GroebnerBasis([M*MAPKK*k1*k2 + Mp_MKP_*(-h5*k2 - h5*k_1), M*MKP*h_6 - M_MKP*h6 + Mp_MKP_*h5, M*Mp_MKP_*(h5*h_6 + h_4*h_6) - M_MKP*Mp*h4*h6 + Mp*Mp_MKP_*h4*h5, M*Mpp_MKP*(h2*k1*k2*k4 + h2*k1*k2*k_3) + Mp*Mp_MKP_*(-h5*k2*k3*k4 - h5*k3*k4*k_1), MAPKK*M_MKP*(h5*h6*k1*k2*k3*k4 + h6*h_4*k1*k2*k3*k4) + MKP*Mp_MKP_*(-h5**2*h_6*k2*k3*k4 - h5**2*h_6*k3*k4*k_1 - h5*h_4*h_6*k2*k3*k4 - h5*h_4*h_6*k3*k4*k_1) + MKP*Mpp_MKP*(-h2*h4*h5*k1*k2*k4 - h2*h4*h5*k1*k2*k_3), MAPKK*Mp*k3*k4 + Mpp_MKP*(-h2*k4 - h2*k_3), MAPKK*Mp_MKP_*(h5*k3*k4 + h_4*k3*k4) + MKP*Mpp_MKP*(-h2*h4*k4 - h2*h4*k_3), MKP*Mp*h4 + Mp_MKP_*(-h5 - h_4), MKP*Mpp*h1 + Mpp_MKP*(-h2 - h_1), M_MAPKK*k2 - Mp_MKP_*h5, M_MKP*Mpp*(h1*h2*h6*k1*k2*k4 + h1*h2*h6*k1*k2*k_3) + Mp*Mp_MKP_*(-h2*h5*h_6*k2*k3*k4 - h2*h5*h_6*k3*k4*k_1 - h5*h_1*h_6*k2*k3*k4 - h5*h_1*h_6*k3*k4*k_1) + Mp_MKP_*Mpp*(-h1*h2*h5*k1*k2*k4 - h1*h2*h5*k1*k2*k_3), M_MKP*Mpp_MKP*(h2*h4*h6*k1*k2*k4 + h2*h4*h6*k1*k2*k_3) + Mp_MKP_**2*(-h5**2*h_6*k2*k3*k4 - h5**2*h_6*k3*k4*k_1 - h5*h_4*h_6*k2*k3*k4 - h5*h_4*h_6*k3*k4*k_1) + Mp_MKP_*Mpp_MKP*(-h2*h4*h5*k1*k2*k4 - h2*h4*h5*k1*k2*k_3), Mp*Mpp_MKP*(h2*h4 + h4*h_1) + Mp_MKP_*Mpp*(-h1*h5 - h1*h_4), Mp_MAPKK*k4 - Mpp_MKP*h2, Mp_MKP*h3*h4 + Mp_MKP_*(-h5*h_3 - h_3*h_4) - Mpp_MKP*h2*h4*uVol], M, MAPKK, MKP, M_MAPKK, M_MKP, Mp, Mp_MAPKK, Mp_MKP, Mp_MKP_, Mpp, Mpp_MKP, domain='ZZ[h1,h2,h3,h4,h5,h6,k1,k2,k3,k4,h_1,h_3,h_4,h_6,k_1,k_3,uVol]', order='lex')

Deficiency, reversibility and linkage classes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can check if the network is weakly reversible:

.. code:: python

    >>> crn.is_weakly_rev
    False

Other features provided by the CRN class are the calculation of the network deficiency,
linkage classes, and terminal complexes
(the following is example S7 in [9]_):

.. code:: python

    >>> net = from_react_strings(["X <-> A", "A -> Ap", "Ap <-> Xp",
    ...                           "Xp + Y <-> B", "B -> Bp", "Bp <-> X + Yp",
    ...                           "Yp + A <-> C", "C -> Cp", "Cp <-> A + Y"])
    >>> net.deficiency
    1
    >>> net.strong_linkage_classes
    [[X, A], [Ap, Xp], [Xp + Y, B], [Bp, X + Yp], [A + Yp, C], [Cp, A + Y]]
    >>> net.linkage_classes
    [[X, A, Ap, Xp], [Xp + Y, B, Bp, X + Yp], [A + Yp, C, Cp, A + Y]]
    >>> net.terminal_complexes
    [Ap, Xp, Bp, X + Yp, Cp, A + Y]
    >>> net.non_terminal_complexes
    [X, A, Xp + Y, B, A + Yp, C]

Other features
~~~~~~~~~~~~~~

*acr_species* looks for species that exhibit absolute concentration robustness using the algorithm in [9]_:

.. code:: python

    >>> net.acr_species()
    ['Yp']

The same method used with the option *subnets = True* will attempt to find a decomposition of the network
in subnetworks, using the network elementary modes, and to use this decomposition to
find species with absolute concentration robustness. Consider example S30 in [9]_:

.. code:: python

    >>> net = from_react_strings(["A + B -> 2B", "B -> A", "2A <-> C", "A + C <-> D"])
    >>> net.acr_species()
    ['A']
    >>> net.acr_species(subnets = True)
    ['A', 'C', 'D']

The influence matrix and adjacency matrix for the directed species reaction graph (as defined in [4]_)
can be created with *influence_matrix()* and *dsr_graph_adj()* respectively:

.. code:: python

    >>> crn = from_react_file("examples/data/reactions/dsr-graph/pos_loops_main")
    >>> crn.influence_matrix(var = "a")
    Matrix([
    [a1_1,    0,    0,    0],
    [   0, a2_2, a2_3,    0],
    [   0,    0, a3_3, a3_4]])
    >>> crn.print_influence_matrix()
            r1     r2     r3     r4
    x1 | g_1_1      0      0      0 |
    x2 |     0  g_2_2  g_2_3      0 |
    x3 |     0      0  g_3_3  g_3_4 |
    >>> crn.dsr_graph_adj()
    Matrix([
    [ 0,  0,  0, 1, 0, 0, 0],
    [ 0,  0,  0, 0, 1, 1, 0],
    [ 0,  0,  0, 0, 0, 1, 1],
    [-1,  1,  0, 0, 0, 0, 0],
    [ 1, -1,  0, 0, 0, 0, 0],
    [ 0, -1,  1, 0, 0, 0, 0],
    [ 0,  1, -1, 0, 0, 0, 0]])

Reduction
---------

The tool offers some methods for the structural reduction of chemical reaction network
and the derivation of kinetic rates.

In the following example, we consider the one-substrate enzyme reaction mechanism,
and eliminate the intermediate *ES* using quasi-steady state approximation ([2]_, [6]_, [8]_):

.. code:: python

    >>> crn = from_sbml("examples/data/sbml/enzyme.xml")
    >>> crn.reactions
    (veq: E + S ->(_comp*veq_kon) ES, veq_rev: ES ->(_comp*veq_koff) E + S, vcat: ES ->(_comp*vcat_kcat) E + P)
    >>> crn.qss('ES')
    >>> for r in crn.reactions: print(r)
    ... 
    veq_vcat: E + S ->(comp*vcat_kcat*veq_kon/(vcat_kcat + veq_koff)) E + P

We can now use a conservation to eliminate the enzyme, and check the new dynamics:

.. code:: python

    >>> from conslaw import ConsLaw
    >>> crn.remove_by_cons('E', ConsLaw('E + ES', 'Et'))
    >>> for r in crn.reactions: print(r)
    ... 
    veq_vcat: S ->(comp*et*vcat_kcat*veq_kon/(s*veq_kon + vcat_kcat + veq_koff)) p
    >>> crn.print_equations()
    dP/dt = comp*Et*S*vcat_kcat*veq_kon/(S*veq_kon + vcat_kcat + veq_koff)
    dS/dt = -comp*Et*S*vcat_kcat*veq_kon/(S*veq_kon + vcat_kcat + veq_koff)

In alternative, we could eliminate the constant species:

.. code:: python

    >>> crn = from_sbml("examples/data/sbml/enzyme.xml")
    >>> crn.qss('ES')
    >>> crn.constant_species
    ['E']
    >>> crn.remove_all_constants()
    >>> for r in crn.reactions: print(r)
    ... 
    veq_vcat: S ->(comp*E*vcat_kcat*veq_kon/(vcat_kcat + veq_koff)) P

or use a rapid equilibrium approximation ([2]_, [6]_, [8]_):

.. code:: python

    >>> crn = from_sbml("examples/data/sbml/enzyme.xml")
    >>> crn.rapid_eq(('ES', 'E + S'), cons_law = ('E', ConsLaw('E + ES', 'Et')))
    >>> for r in crn.reactions: print(r)
    ... 
    vcat: S ->(comp*Et*vcat_kcat*veq_kon/(S*veq_kon + veq_koff)) P

With the method *remove* we can use a combination of the reduction methods:

.. code:: python

    >>> bi_uni_random.remove(rapid_eq = [('ea', 'e + a'), ('eb', 'e + b')], 
                             qss = ['eab'], 
                             cons_law = ('e', ConsLaw('e + ea + eb + eab', 'et')))
    >>> for r in bi_uni_random.reactions: print(r)
    ... 
    r2_r4: a + b ->(et*k1*k3*k5*k_2/(a*b*k1*k3*k_2 + a*b*k2*k4*k_1 + a*k1*k5*k_2 + a*k1*k_2*k_3 + a*k1*k_2*k_4 + b*k2*k5*k_1 + b*k2*k_1*k_3 + b*k2*k_1*k_4 + k5*k_1*k_2 + k_1*k_2*k_3 + k_1*k_2*k_4)) p
    r3_r4: a + b ->(et*k2*k4*k5*k_1/(a*b*k1*k3*k_2 + a*b*k2*k4*k_1 + a*k1*k5*k_2 + a*k1*k_2*k_3 + a*k1*k_2*k_4 + b*k2*k5*k_1 + b*k2*k_1*k_3 + b*k2*k_1*k_4 + k5*k_1*k_2 + k_1*k_2*k_3 + k_1*k_2*k_4)) p

We can merge reactions with the same reactant and product:

.. code:: python

    >>> bi_uni_random.merge_reactions()
    >>> for r in bi_uni_random.reactions: print(r)
    ... 
    r2_r4r3_r4: a + b ->(et*k5*(k1*k3*k_2 + k2*k4*k_1)/(a*b*k1*k3*k_2 + a*b*k2*k4*k_1 + a*k1*k5*k_2 + a*k1*k_2*k_3 + a*k1*k_2*k_4 + b*k2*k5*k_1 + b*k2*k_1*k_3 + b*k2*k_1*k_4 + k5*k_1*k_2 + k_1*k_2*k_3 + k_1*k_2*k_4)) p

Saving models
-------------

Chemical reaction networks can be saved to SBML files

.. code:: python

    >>> crn.save_sbml("examples/data/sbml/enzyme_simplified.xml")

or as reaction files (by default the strings representing reactions contain the kinetic parameters;
use *rate = True* to save the reaction rates instead):

.. code:: python

    >>> crn.save_reaction_file("examples/data/reactions/enzyme_simplified", rate = True)

References
----------

.. [1] Angeli, D. (2009). *A tutorial on Chemical Reaction Networks dynamics*. In Control Conference (ECC), 2009 European (pp. 649-657). IEEE.

.. [2] Cornish-Bowden, A. (1987). *Fundamentals of Enzyme Kinetics*. Elsevier Science.

.. [3] Feinberg, M. (1979). *Lectures on chemical reaction networks*. Notes of lectures given at the Mathematics Research Center, University of Wisconsin.

.. [4] Feliu, E., & Wiuf, C. (2015). *Finding the positive feedback loops underlying multi-stationarity*. BMC systems biology, 9(1), 1.

.. [5] Gunawardena, J. (2003). *Chemical reaction network theory for in-silico biologists*, http://vcp.med.harvard.edu/papers/crnt.pdf.

.. [6] Ingalls, Brian. (2013). *Mathematical Modelling in Systems Biology: An Introduction.*, https://www.math.uwaterloo.ca/~bingalls/MMSB/.

.. [7] Juty, N., et al. (2015). *BioModels: content, features, functionality, and use.* CPT: pharmacometrics \& systems pharmacology, 4(2), pp.55-68.

.. [8] Segel, I. H. (1975). *Enzyme kinetics*. Vol. 957. Wiley, New York.

.. [9] Shinar, G., Feinberg, M. (2010), *Structural sources of robustness in biochemical reaction networks*, Science.
