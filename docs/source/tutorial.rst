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
stoichiometric matric *stoich\_matrix*, the matrix of stoichiometric
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


Check if the network is weakly reversible:

.. code:: python

    >>> crn.is_weakly_rev
    False

Other features provided by the CRN class are the calculation of the network deficiency,
linkage classes, and terminal complexes
(the following is example S7 in [1]_):

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

*acr_species* looks for species that exhibit absolute concentration robustness using the algorithm in [1]_:

.. code:: python

    >>> net.acr_species()
    ['Yp']

The same method used with the option *subnets = True* will attempt to find a decomposition of the network
in subnetworks, using the network elementary modes, and to use this decomposition to
find species with absolute concentration robustness. For example S30 in [1]_:

.. code:: python

    >>> net = from_react_strings(["A + B -> 2B", "B -> A", "2A <-> C", "A + C <-> D"])
    >>> net.acr_species()
    ['A']
    >>> net.acr_species(subnets = True)
    ['A', 'C', 'D']


Reduction
~~~~~~~~~

Eliminate the intermediate *ES* using quasi-steady state approximation:

.. code:: python

    >>> crn = from_sbml("examples/data/sbml/enzyme.xml")
    >>> crn.qss('ES')
    >>> for r in crn.reactions: print(r)
    ... 
    veq_vcat: E + S ->(comp*vcat_kcat*veq_kon/(vcat_kcat + veq_koff)) E + P

Use a conservation to eliminate the enzyme, and check the new dynamics:

.. code:: python

    >>> from conslaw import ConsLaw
    >>> crn.remove_by_cons('E', ConsLaw('E + ES', 'Et'))
    >>> for r in crn.reactions: print(r)
    ... 
    veq_vcat: S ->(comp*et*vcat_kcat*veq_kon/(s*veq_kon + vcat_kcat + veq_koff)) p
    >>> crn.print_equations()
    dP/dt = comp*Et*S*vcat_kcat*veq_kon/(S*veq_kon + vcat_kcat + veq_koff)
    dS/dt = -comp*Et*S*vcat_kcat*veq_kon/(S*veq_kon + vcat_kcat + veq_koff)

In alternative, eliminate the constant species:

.. code:: python

    >>> crn = from_sbml("examples/data/sbml/enzyme.xml")
    >>> crn.qss('ES')
    >>> crn.constant_species
    ['e']
    >>> crn.remove_all_constants()
    >>> for r in crn.reactions: print(r)
    ... 
    veq_vcat: S ->(comp*E*vcat_kcat*veq_kon/(vcat_kcat + veq_koff)) P

Use rapid equilibrium instead (and the conservation law):

.. code:: python

    >>> crn = from_sbml("examples/data/sbml/enzyme.xml")
    >>> crn.rapid_eq(('ES', 'E + S'), cons_law = ('E', ConsLaw('E + ES', 'Et')))
    >>> for r in crn.reactions: print(r)
    ... 
    vcat: S ->(comp*Et*vcat_kcat*veq_kon/(S*veq_kon + veq_koff)) P

Use a combination of the reduction methods:

.. code:: python

    >>> bi_uni_random.remove(rapid_eq = [('ea', 'e + a'), ('eb', 'e + b')], 
                           qss = ['eab'], 
                           cons_law = ('e', ConsLaw('e + ea + eb + eab', 'et')))
    >>> for r in bi_uni_random.reactions: print(r)
    ... 
    r2_r4: a + b ->(et*k1*k3*k5*k_2/(a*b*k1*k3*k_2 + a*b*k2*k4*k_1 + a*k1*k5*k_2 + a*k1*k_2*k_3 + a*k1*k_2*k_4 + b*k2*k5*k_1 + b*k2*k_1*k_3 + b*k2*k_1*k_4 + k5*k_1*k_2 + k_1*k_2*k_3 + k_1*k_2*k_4)) p
    r3_r4: a + b ->(et*k2*k4*k5*k_1/(a*b*k1*k3*k_2 + a*b*k2*k4*k_1 + a*k1*k5*k_2 + a*k1*k_2*k_3 + a*k1*k_2*k_4 + b*k2*k5*k_1 + b*k2*k_1*k_3 + b*k2*k_1*k_4 + k5*k_1*k_2 + k_1*k_2*k_3 + k_1*k_2*k_4)) p

Merge reactions with the same reactant and product:

.. code:: python

    >>> bi_uni_random.merge_reactions()
    >>> for r in bi_uni_random.reactions: print(r)
    ... 
    r2_r4r3_r4: a + b ->(et*k5*(k1*k3*k_2 + k2*k4*k_1)/(a*b*k1*k3*k_2 + a*b*k2*k4*k_1 + a*k1*k5*k_2 + a*k1*k_2*k_3 + a*k1*k_2*k_4 + b*k2*k5*k_1 + b*k2*k_1*k_3 + b*k2*k_1*k_4 + k5*k_1*k_2 + k_1*k_2*k_3 + k_1*k_2*k_4)) p

Saving models
~~~~~~~~~~~~~

Chemical reaction networks can be saved to SBML files, and reaction files:

.. code:: python

    >>> crn.save_sbml("examples/data/sbml/enzyme_simplified.xml")
    >>> crn.save_reaction_file("examples/data/reactions/enzyme_simplified")

Other features
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

Groebner basis...

References
==========

.. [1] Shinar, G., Feinberg, M. (2010), *Structural sources of robustness in biochemical reaction networks*, Science.
