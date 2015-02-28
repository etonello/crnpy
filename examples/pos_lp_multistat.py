#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

from functools import reduce
from itertools import combinations
import networkx as nx
from operator import mul

from crnpy.crn import from_react_file

# Implementing algorithm for finding positive
# feedback loops generating multistationarity, from
# Feliu, E., & Wiuf, C. (2015).
# Finding the positive feedback loops underlying multi-stationarity.
# BMC systems biology, 9(1), 1


def cycle_sign(c, graph):
    return reduce(mul, [graph.get_edge_data(*d)['weight'] for d in zip(c, c[1:] + [c[0]])], 1)


def nucleus_sign(nucleus, graph):
    return (-1)**sum([1 for c in nucleus if (len(c)/2)%2==0])*reduce(mul, [cycle_sign(c, graph) for c in nucleus], 1)


def pos_fb_loops_multistat(crn, debug = False):
    A = crn.stoich_matrix
    im = crn.influence_matrix()
    M = A * im.T
    # Mtilde: replace the i_j row of M with j-th element of the basis
    omega = crn._cons_laws()
    for ij in omega:
        M[list(ij).index(1), :] = ij.T
    p = M.det()
    Gfull = nx.from_numpy_matrix(crn.dsr_graph_adj(), create_using = nx.DiGraph())

    if debug:
        for r in crn.reactions: print(r)
        print("")
        print("Stoichiometric matrix")
        crn.print_stoich_matrix()
        print("")
        print("Conservations: {}".format(crn.cons_laws))
        print("")
        print("Influence matrix")
        print_matrix(im, crn.species, [r.reactionid for r in crn.reactions])
        print("")
        print("M = A Z^t")
        print_matrix(M, crn.species, crn.species)
        print("")
        print("Determinant:")
        print(p)
        print("")

    if p == 0:
        print("Determinant is zero.")
        return []

    # if both not zero, there are multiple signs
    if p + p.coeff(-1) != 0 and p.coeff(-1) != 0:
        print("Multiple signs: cannot exclude multistationarity")
    else:
        print("Same sign: multistationarity excluded")
        return []

    s = A.rank()
    sign = (-1)**(s + 1)
    if debug: print("Rank: {}, sign: {}".format(s, sign))
    # select the terms with sign sign
    sumterms = p.coeff(-1) if sign == -1 else p + p.coeff(-1)

    if sumterms.func.__name__ == 'Add':
        terms = sumterms.args
    else:
        # if there is only one addend
        terms = [sumterms]

    nodes = crn.species + crn.reactionids
    # find positive feedback loops associated to terms with sign s
    pfls = []
    if debug: print("Looking for {} nuclei".format(2 * s))
    for term in terms:
        if debug: print("Term: {}".format(term))

        # create graph with elements of Z only involved in term
        G = nx.from_numpy_matrix(crn.dsr_graph_adj(keep = term.args), create_using = nx.DiGraph())

        # find the elementary cycles
        cycles = list(nx.simple_cycles(G))

        if debug:
            print("Elementary cycles:")
            for c in cycles: print([nodes[i] for i in c], cycle_sign(c, G))

        # look at all cycle combinations,
        # find those with the right length (2s),
        # check if the sign is s;
        # if so, find the positive feedback loops within
        for t in range(1, s+1):
            for c in list(combinations(cycles, t)):
                if sum([len(d) for d in c]) == 2 * s:
                    if debug: print(2*s, "nucleus", [[nodes[i] for i in d] for d in c], nucleus_sign(c, Gfull))
                    csign = nucleus_sign(c, Gfull)
                    if csign == sign:
                        for d in c:
                            if cycle_sign(d, Gfull) == 1:
                                pfl = [nodes[i] for i in d]
                                if pfl not in pfls: pfls.append(pfl)
    return pfls


def pfloops_examples(debug = False):
    folder = "data/reactions/dsr-graph/"
    filenames = ["pos_loops_main",
                 "pos_loops_main_v2",
                 "ubiquitination",
                 "phosphorylation",
                 "signalling_cascades",
                 "apoptosis"]

    for i in range(len(filenames)):
        crn = from_react_file(os.path.join(folder, filenames[i]))
        print(filenames[i])
        for r in crn.reactions: print(r)
        pfls = pos_fb_loops_multistat(crn, debug)
        if len(pfls) > 0:
            print("Positive feedback loops found:")
            for pfl in pfls: print(pfl)
        print("")


if __name__ == "__main__":
    pfloops_examples()
