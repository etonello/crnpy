import sys
import os
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

from crnpy.crn import CRN, from_react_file

# G-protein signalling
# Ingalls, Brian. "Mathematical Modelling in Systems Biology: An Introduction.", 2013.
# 6.1.2

print("Creating model...")
crn = from_react_file("data/reactions/g-protein_signalling")
crn.inspect(True)

print("Assuming ligand is constant:")
crn.remove_constant('l')
crn.inspect(True)
