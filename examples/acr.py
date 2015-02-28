import sys
import os
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

from crnpy.crn import CRN, from_react_file


filename = "data/reactions/acr/acr_1"
crn = from_react_file(filename)
for r in crn.reactions: print(r)
print
print("Terminal complexes: {}".format(crn.terminal_complexes()))
print("Nonterminal complexes: {}".format(crn.non_terminal_complexes()))
print("ACR in species: {}".format(", ".join(crn.acr_species())))

print
# print latex
for r in crn.reactions:
    print(r.latex() + " \\\\")
print
