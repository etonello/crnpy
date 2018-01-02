import sys
import os
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

from crnpy.crn import CRN, from_react_file

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"

# Cooperative binding

print("Creating model...")
crn = from_react_file("data/reactions/cooperative_binding")
crn.inspect(True)

print("")

print("Remove ps1, ps2 and ps3 by qssa")
crn.remove(qss = ['ps1', 'ps2', 'ps3'])
for s, f in crn.removed_species: print(s + " = " + str(f))
crn.inspect(True)
