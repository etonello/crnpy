#!/usr/bin/python
# -*- coding: latin-1 -*-

import sys
import os
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

from crnpy.crn import *

from itertools import product
from pprint import pprint
from pulp import *


# Algorithm for computing sparse and dense realizations of CRNs based on
# "Computing sparse and dense realizations of reaction kinetic systems", Gábor Szederkényi (2010)
# Journal of Mathematical Chemistry, 47(2), 551-568.


def sparse_dense_realiz(crn, params, debug = False, sparse = True, lp_model_file = None):
    m = crn.n_complexes
    n = crn.n_species
    Y = crn.complex_matrix
    Ak = -crn.laplacian
    M = Y*Ak

    prob = LpProblem("Realization", LpMinimize)

    # upper bound for elements of A
    u = 1e+3


    ### Variables ###
    offdiagonal = [(i, j) for i, j in product(range(m), repeat = 2) if i != j]

    # Kirkoff matrix to identify
    A = LpVariable.dicts('A', offdiagonal, 0, u)

    # Variables to count number of reactions
    delta = LpVariable.dicts('delta', offdiagonal, 0, 1, 'Integer')


    ### Objective ###
    prob += lpSum(delta) * (1 if sparse else -1)


    ## Constraints ###

    # Y*A = M
    for s in range(n):
        for c in range(m):
            prob += lpSum((Y[s, k] - Y[s, c])*A[(k, c)] for k in range(m) if k != c) == M[s, c]

    # delta counts number of reactions
    for p in offdiagonal:
        prob += delta[p] >= A[p] * 1./ (u + 1), "Count reaction " + str(p)
        prob += delta[p] <= A[p] * u, "No reaction " + str(p)


    ### Solve ###
    if lp_model_file:
        prob.writeLP("test.lp")
    msg = 0
    if debug: msg = 1
    status = prob.solve(GLPK(msg = msg))
    if status != 1:
        print("Unable to solve. Status: {}".format(status))


    ## Solution
    sol = dict((v.name, v.varValue) for v in prob.variables() if v.varValue != None)
    if debug:
        print("Objective = {}".format(value(prob.objective)))
        print("Solution: {}".format(sol))

    reacts = []
    for c1, c2 in offdiagonal:
        param = sol["A_({},_{})".format(c1, c2)]
        if param != 0:
            reacts.append(Reaction("r_{}_{}".format(c2, c1), # reactionid
                                   crn.complexes[c2], # reactant
                                   crn.complexes[c1], # product
                                   param*crn.complexes[c2].ma())) # rate

    crn_sol = from_reacts(reacts)

    return crn_sol


def compare_realizations(crn, params):
    # find sparse and dense realizations for the given parameters
    print("Original network: {} reactions".format(crn.n_reactions))
    crn.print_laplacian()
    for e in crn.format_equations(): print (e)
    print(crn.format_deficiency())
    print("Weakly reversible: {}".format(crn.is_weakly_rev))
    print('')

    # replace the kinetic parameters with params
    crn.set_params(params)

    crn.print_laplacian(numeric = True)
    pprint(crn.reactions)
    for e in crn.format_equations(): print (e)
    print('')

    sparse = sparse_dense_realiz(crn, params, debug = False, sparse = True)
    print("Sparse realization: {} reactions".format(sparse.n_reactions))
    sparse.print_laplacian(numeric = True, precision = 4)
    pprint(sparse.reactions)
    for e in sparse.format_equations(): print (e)
    print(sparse.format_deficiency())
    print("Weakly reversible: {}".format(sparse.is_weakly_rev))
    print('')

    dense = sparse_dense_realiz(crn, params, debug = False, sparse = False)
    print("Dense realization: {} reactions".format(dense.n_reactions))
    dense.print_laplacian(numeric = True, precision = 4)
    pprint(dense.reactions)
    for e in dense.format_equations(): print (e)
    print(dense.format_deficiency())
    print("Weakly reversible: {}".format(dense.is_weakly_rev))


if __name__ == "__main__":
    # Example 2.1
    print("Example 2.1")
    crn = from_react_file("data/reactions/realizations/ex_2.1_szed_2010")
    # values of the kinetic parameters
    params = {'k1': 1,
              'k2': 1.1,
              'k3': 1,
              'k4': 1,
              'k5': 1.1,
              'k6': 0.1,
              'k7': 3,
              'k8': 1}
    compare_realizations(crn, params)

    print('')

    # Example 4.2
    print("Example 4.2")
    crn = from_react_file("data/reactions/realizations/ex_4.2_szed_2010")
    # values of the kinetic parameters
    params = dict((k, 1) for k in crn.kinetic_params)
    compare_realizations(crn, params)

    print('')

    # Example 4.3
    print("Example 4.3")
    crn = from_react_file("data/reactions/realizations/ex_4.3_szed_2010")
    # values of the kinetic parameters
    params = dict((k, 1) for k in crn.kinetic_params)
    compare_realizations(crn, params)
