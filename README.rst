crnpy
=====

crnpy is a python library for the manipulation and analysis of chemical
reaction networks.

Install
-------

crnpy can be installed from source. After downloading the repository, run

      $ python setup.py install --user

crnpy requires
`libSBML <http://sbml.org/Software/libSBML>`_,
`SciPy <https://www.scipy.org/scipylib/index.html>`_,
`NumPy <http://www.numpy.org/>`_ and
`SymPy <http://www.sympy.org/en/index.html>`_.
To generate the network invariants, `pycddlib <http://pycddlib.readthedocs.io/en/latest/quickstart.html#installation>`_ is also required.

Some example scripts require `NetworkX <http://networkx.github.io/>`_ or `PuLP <http://pythonhosted.org/PuLP/>`_.

Getting Started
---------------

We can create a network from an SBML file

.. code:: python

      >>> from crnpy.crn import CRN, from_sbml, from_react_strings, from_react_file
      >>> crn = from_sbml("examples/data/sbml/enzyme.xml")

from a file containing a list of reactions in human-readable format

.. code:: python

      >>> crn = from_react_file("examples/data/reactions/biomodels/biomd0000000026")

or directly from a list of reaction strings:

.. code:: python

      >>> crn = from_react_strings(["A <-> B", "2A + C <-> D", "D -> E", "E -> 2A + C"])

Now we can explore some properties of the network. For example, we can
look at the stoichiometric matrix

.. code:: python

      >>> crn.print_stoich_matrix()
          r0  r0_rev  r1  r1_rev  r2  r3
      A | -1       1  -2       2   0   2 |
      B |  1      -1   0       0   0   0 |
      C |  0       0  -1       1   0   1 |
      D |  0       0   1      -1  -1   0 |
      E |  0       0   0       0   1  -1 |

at the derivatives of the species concentrations and the conservation laws:

.. code:: python

      >>> crn.stoich_matrix * crn.rates
      Matrix([
      [-2*A**2*C*k_r1 - A*k_r0 + B*k_r0_rev + 2*D*k_r1_rev + 2*E*k_r3],
      [                                           A*k_r0 - B*k_r0_rev],
      [                            -A**2*C*k_r1 + D*k_r1_rev + E*k_r3],
      [                             A**2*C*k_r1 - D*k_r1_rev - D*k_r2],
      [                                               D*k_r2 - E*k_r3]])
      >>> crn.cons_laws
      (A + B + 2*D + 2*E, C + D + E)

We can check whether the conditions of the deficiency zero theorem are
satisfied:

.. code:: python

      >>> crn.is_ma
      True
      >>> crn.deficiency
      0
      >>> crn.is_weakly_rev
      True

For more information, a tutorial is available as well as some
example scripts.


Citation
--------

If you use crnpy for your work, please cite

Elisa Tonello, CrnPy: a python library for the analysis of chemical reaction networks, 2016.
