import sys
import os
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

from crnpy.crn import CRN, from_react_file
from crnpy.conslaw import ConsLaw

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"

# G-protein signalling
# Ingalls, Brian. "Mathematical Modelling in Systems Biology: An Introduction.", 2013.
# 6.1.2

print("Creating model...")
crn = from_react_file("data/reactions/g-protein_signalling")
crn.inspect(True, invariants = True)

print("")

print("Assuming ligand is constant:")
crn.remove_constant('L')
crn.inspect(True, invariants = True)

print("")

print("Using conservations to reduce the number of equations:")
crn.remove_by_cons('RL', ConsLaw('R + RL', 'Rtot'))
crn.remove_by_cons('Gbg', ConsLaw('G + Gbg', 'Gtot1'))
crn.remove_by_cons('Gd', ConsLaw('G + Ga + Gd', 'Gtot'))
crn.inspect(True, invariants = True)
