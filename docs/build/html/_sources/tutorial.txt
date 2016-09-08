crnpy tutorial
==============

Defining a chemical reaction network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A chemical reaction network (CRN) object can be created from an SBML
file:

.. code:: python

    >>> from crnpy.crn import CRN, from_sbml, from_react_strings, from_react_file
    >>> crn = from_sbml("examples/data/sbml/enzyme.xml")

or a list of reactions:

.. code:: python

    >>> reactions = ['e + a (k_1)<->(k1) ea',
    ...              'e + b (k_2)<->(k2) eb',
    ...              'ea + b (k_3)<->(k3) eab',
    ...              'eb + a (k_4)<->(k4) eab',
    ...              'eab ->(k5) e + p']
    >>> bi_uni_random = from_react_strings(reactions)

The list of reactions can be read from a text file, using
*from\_react\_file*:

.. code:: python

    >>> another_net = from_react_file(reactionfilename)

A string specifying a reaction must contain either the characters "->"
for unidirectional reactions, or the characters "<->" for reversible
reactions. Reversible reactions will be converted to two separate
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
by a colon. For example, we can specify the last reaction as

.. code:: python

    >>> one_react_net = from_react_strings(["r_in: 2a + b ->(k*a**2*b/(a+b+c)) 3 c"], rate = True)

If an id is not specified, reactions are assigned an id of the form
*r\_n*, with n an integer, starting from 0. If a reversible reaction is
defined, for example

.. code:: python

    >>> rev_react_net = from_react_strings(["r1: a + b <-> c"])

then two reactions will be created, one with id r1 with a + b as a
reactant and c as product, and one with id r1\_rev with c as reactant
and a + b as product. As shown in the last example, kinetic parameters
are optional. In the same example, reaction r1 is assigned a parameter
symbol k\_r1, and the reverse reaction is assigned the parameter symbol
k\_r1\_rev:

.. code:: python

    >>> rev_react_net.reactions
    (r1: a + b ->(k_r1) c, r1_rev: c ->(k_r1_rev) a + b)

Comments can be added to a reaction file using the symbol "#". Anything
appearing after the hash sign will ignored.

Exploring chemical reaction networks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The library implements some elements of chemical reaction network theory ([1]_, [3]_, [4]_).

Attributes of a CRN object include the network species, complexes, and
reactions:

.. code:: python

    >>> crn.species, crn.complexes
    (('E', 'ES', 'P', 'S'), (E + S, ES, E + P))

    >>> for r in crn.reactions: print(r)
    ... 
    veq: E + S ->(comp*veq_kon) ES
    veq_rev: ES ->(comp*veq_koff) E + S
    vcat: ES ->(comp*vcat_kcat) E + P

    >>> for r in bi_uni_random.reactions: print(r)
    ... 
    r0: a + e ->(k1) ea
    r0_rev: ea ->(k_1) a + e
    r1: b + e ->(k2) eb
    r1_rev: eb ->(k_2) b + e
    r2: b + ea ->(k3) eab
    r2_rev: eab ->(k_3) b + ea
    r3: a + eb ->(k4) eab
    r3_rev: eab ->(k_4) a + eb
    r4: eab ->(k5) e + p

Available matrices associated to the reaction network are the
stoichiometric matrix *stoich\_matrix*, the matrix of stoichiometric
coefficients *complex\_matrix* (often called Y in the literature), the
Laplacian of the graph of complexes *laplacian*, and its negation *kinetic_matrix*,
the incidence matrix of the complex graph *incidence\_matrix*.

.. code:: python

    >>> crn.stoich_matrix()
    Matrix([
    [-1,  1,  1],
    [ 1, -1, -1],
    [ 0,  0,  1],
    [-1,  1,  0]])

Special methods are available to print some matrices. For example, for
the stoichiometry matrix:

.. code:: python

      >>> bi_uni_random.print_stoich_matrix()
          r0  r0_rev  r1  r1_rev  r2  r2_rev  r3  r3_rev  r4
    a   | -1       1   0       0   0       0  -1       1   0 |
    b   |  0       0  -1       1  -1       1   0       0   0 |
    e   | -1       1  -1       1   0       0   0       0   1 |
    ea  |  1      -1   0       0  -1       1   0       0   0 |
    eab |  0       0   0       0   1      -1   1      -1  -1 |
    eb  |  0       0   1      -1   0       0  -1       1   0 |
    p   |  0       0   0       0   0       0   0       0   1 |

We can look for example at the system of ODEs associated to the network,
and at the conservation laws:

.. code:: python

    >>> crn.print_equations()
    dE/dt = -comp*E*S*veq_kon + comp*ES*vcat_kcat + comp*ES*veq_koff
    dES/dt = comp*E*S*veq_kon - comp*ES*vcat_kcat - comp*ES*veq_koff
    dP/dt = comp*ES*vcat_kcat
    dS/dt = -comp*E*S*veq_kon + comp*ES*veq_koff
    >>> crn.cons_laws
    (E - P - S, ES + P + S)

or get a list of intermediate species:

.. code:: python

    >>> crn.intermediate_species
    ['E', 'ES', 'S']


We can check if the network is weakly reversible:

.. code:: python

    >>> crn.is_weakly_rev
    False

Other features provided by the CRN class are the calculation of the network deficiency,
linkage classes, and terminal complexes
(the following is example S7 in [8]_):

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
    >>> net.deficiency
    1
    >>> net.terminal_complexes
    [Ap, Xp, Bp, X + Yp, Cp, A + Y]
    >>> net.non_terminal_complexes
    [X, A, Xp + Y, B, A + Yp, C]

*acr_species* looks for species that exhibit absolute concentration robustness using the algorithm in [8]_:

.. code:: python

    >>> net.acr_species()
    ['Yp']

The same method used with the option *subnets = True* will attempt to find a decomposition of the network
in subnetworks, using the network elementary modes, and to use this decomposition to
find species with absolute concentration robustness. Consider example S30 in [8]_:

.. code:: python

    >>> net = from_react_strings(["A + B -> 2B", "B -> A", "2A <-> C", "A + C <-> D"])
    >>> net.acr_species()
    ['A']
    >>> net.acr_species(subnets = True)
    ['A', 'C', 'D']


Reduction
~~~~~~~~~

The tool offers some methods for the structural reduction of chemical reaction network
and the derivation of kinetic rates (the algorithms used in the following examples are described in [9]_).

In the following example, we consider the one-substrate enzyme reaction mechanism,
and eliminate the intermediate *ES* using quasi-steady state approximation ([2]_, [5]_, [7]_):

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

or use a rapid equilibrium approximation ([2]_, [5]_, [7]_):

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
~~~~~~~~~~~~~

Chemical reaction networks can be saved to SBML files

.. code:: python

    >>> crn.save_sbml("examples/data/sbml/enzyme_simplified.xml")

or as reaction files (by default the strings representing reactions contain the kinetic parameters;
use *rate = True* to save the reaction rates instead):

.. code:: python

    >>> crn.save_reaction_file("examples/data/reactions/enzyme_simplified", rate = True)

Other examples
~~~~~~~~~~~~~~

Create a model and look at its deficiency and elementary modes:

.. code:: python

    >>> reactions = ['r1: a ->(k1) b + y',
    ...              'r2: y ->(k2) c',
    ...              'r3: b + c ->(k3) a']
    >>> example = from_react_strings(reactions)
    >>> example.deficiency
    1
    >>> example.n_complexes, example.n_linkage_classes, example.stoich_matrix.rank()
    (5, 2, 2)
    >>> example.elem_modes
    [[1, 1, 0, 1], [0, 1, 1, 0]]
    >>> example.format_elem_modes()
    [[r1 + r2 + r3], [r2 + r2_rev]]
    >>> example.is_weakly_rev
    False

We can check how the elementary modes change if *y* is eliminated:

.. code:: python

    >>> example.qss('y')
    >>> example.reactions
    (r3: b + c ->(k3) a, r1_r2: a ->(k1) b + c)
    >>> example.deficiency
    0
    >>> example.format_elem_modes()
    [[r1_r2 + r3]]
    >>> example.is_weakly_rev
    True

Check if two networks are dynamically equivalent:

.. code:: python

    >>> net1 = from_react_strings(['a ->(k) a + 2b'])
    >>> net2 = from_react_strings(['a ->(2*k) a + b'])
    >>> net1.is_dyn_eq(net2)
    True

We can create a chemical reaction network for a network from the `BioModels <http://biomodels.caltech.edu/>`_ database [6]_:

.. code:: python

    crn = from_react_file("examples/data/reactions/biomodels/biomd0000000026")
    >>> for r in crn.reactions: print(r)
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

and generate a Gröbner basis for the steady state ideal:

.. code:: python

    >>> crn.groebner()
    GroebnerBasis([M*MAPKK*k1*k2 + Mp_MKP_*(-h5*k2 - h5*k_1), M*MKP*h_6 - M_MKP*h6 + Mp_MKP_*h5, M*Mp_MKP_*(h5*h_6 + h_4*h_6) - M_MKP*Mp*h4*h6 + Mp*Mp_MKP_*h4*h5, M*Mpp_MKP*(h2*k1*k2*k4 + h2*k1*k2*k_3) + Mp*Mp_MKP_*(-h5*k2*k3*k4 - h5*k3*k4*k_1), MAPKK*M_MKP*(h5*h6*k1*k2*k3*k4 + h6*h_4*k1*k2*k3*k4) + MKP*Mp_MKP_*(-h5**2*h_6*k2*k3*k4 - h5**2*h_6*k3*k4*k_1 - h5*h_4*h_6*k2*k3*k4 - h5*h_4*h_6*k3*k4*k_1) + MKP*Mpp_MKP*(-h2*h4*h5*k1*k2*k4 - h2*h4*h5*k1*k2*k_3), MAPKK*Mp*k3*k4 + Mpp_MKP*(-h2*k4 - h2*k_3), MAPKK*Mp_MKP_*(h5*k3*k4 + h_4*k3*k4) + MKP*Mpp_MKP*(-h2*h4*k4 - h2*h4*k_3), MKP*Mp*h4 + Mp_MKP_*(-h5 - h_4), MKP*Mpp*h1 + Mpp_MKP*(-h2 - h_1), M_MAPKK*k2 - Mp_MKP_*h5, M_MKP*Mpp*(h1*h2*h6*k1*k2*k4 + h1*h2*h6*k1*k2*k_3) + Mp*Mp_MKP_*(-h2*h5*h_6*k2*k3*k4 - h2*h5*h_6*k3*k4*k_1 - h5*h_1*h_6*k2*k3*k4 - h5*h_1*h_6*k3*k4*k_1) + Mp_MKP_*Mpp*(-h1*h2*h5*k1*k2*k4 - h1*h2*h5*k1*k2*k_3), M_MKP*Mpp_MKP*(h2*h4*h6*k1*k2*k4 + h2*h4*h6*k1*k2*k_3) + Mp_MKP_**2*(-h5**2*h_6*k2*k3*k4 - h5**2*h_6*k3*k4*k_1 - h5*h_4*h_6*k2*k3*k4 - h5*h_4*h_6*k3*k4*k_1) + Mp_MKP_*Mpp_MKP*(-h2*h4*h5*k1*k2*k4 - h2*h4*h5*k1*k2*k_3), Mp*Mpp_MKP*(h2*h4 + h4*h_1) + Mp_MKP_*Mpp*(-h1*h5 - h1*h_4), Mp_MAPKK*k4 - Mpp_MKP*h2, Mp_MKP*h3*h4 + Mp_MKP_*(-h5*h_3 - h_3*h_4) - Mpp_MKP*h2*h4*uVol], M, MAPKK, MKP, M_MAPKK, M_MKP, Mp, Mp_MAPKK, Mp_MKP, Mp_MKP_, Mpp, Mpp_MKP, domain='ZZ[h1,h2,h3,h4,h5,h6,k1,k2,k3,k4,h_1,h_3,h_4,h_6,k_1,k_3,uVol]', order='lex')

References
~~~~~~~~~~

.. [1] Angeli, D. (2009). *A tutorial on Chemical Reaction Networks dynamics*. In Control Conference (ECC), 2009 European (pp. 649-657). IEEE.

.. [2] Cornish-Bowden, A. (1987). *Fundamentals of Enzyme Kinetics*. Elsevier Science.

.. [3] Feinberg, M. (1979). *Lectures on chemical reaction networks*. Notes of lectures given at the Mathematics Research Center, University of Wisconsin.

.. [4] Gunawardena, J. (2003). *Chemical reaction network theory for in-silico biologists*, http://vcp.med.harvard.edu/papers/crnt.pdf.

.. [5] Ingalls, Brian. (2013). *Mathematical Modelling in Systems Biology: An Introduction.*, https://www.math.uwaterloo.ca/~bingalls/MMSB/.

.. [6] Juty, N., et al. (2015). *BioModels: content, features, functionality, and use.* CPT: pharmacometrics \& systems pharmacology, 4(2), pp.55-68.

.. [7] Segel, I. H. (1975). *Enzyme kinetics*. Vol. 957. Wiley, New York.

.. [8] Shinar, G., Feinberg, M. (2010), *Structural sources of robustness in biochemical reaction networks*, Science.

.. [9] Tonello, E., Owen, M. R., Farcot, E. (2016). *On the elimination of intermediate species in chemical reaction networks*.
