#!/usr/bin/env python

"""Derivation of kinetic rates.
References
==========
[1] Cornish-Bowden, Athel. "Fundamentals of enzyme kinetics". Elsevier, 2014.
[2] Ingalls, Brian. "Mathematical Modelling in Systems Biology: An Introduction.", 2013.
[3] SBMLsqueezer documentation, http://www.cogsys.cs.uni-tuebingen.de/software/SBMLsqueezer/doc/KineticLaws2.pdf.
[4] Segel, Irwin H. "Enzyme kinetics". Vol. 957. Wiley, New York, 1975.
"""

import os
import sys
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.realpath(__file__)), '..'))

import sympy as sp
import warnings

from crnpy.conslaw import ConsLaw
from crnpy.crn import CRN, from_sbml, from_react_file, from_react_strings
from crnpy.crncomplex import sympify
from crnpy.parsereaction import parse_reactions, parse_complex

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


input_reactions = os.path.join(os.path.dirname(__file__), 'data', 'reactions')
input_sbml = os.path.join(os.path.dirname(__file__), 'data', 'sbml')


def get_monoms(expr, varlist):
    """Return the monomials in expr with variables in varlist."""
    varlist = list(map(sympify, varlist))
    degrees = sp.Poly(expr, varlist).monoms()
    return [sp.Mul(*(varlist[i]**m[i] for i in range(len(m)))) for m in degrees]


def fail_if_not_equal(a, b):
    if a != b:
        raise ValueError('Rate does not match.')


def eqs_match(origeqs, origspecies, removed, neweqs, newspecies, debug = False):
    """Check that the original equations match the new equations after
    replacing the species with the expression in removed."""
    origeqsdict = dict((origspecies[s], origeqs[s]) for s in range(len(origspecies)))
    replaceeqs = {}
    for s in origeqsdict:
        e = origeqsdict[s]
        for i in range(len(removed)): e = e.subs(removed[i][0], removed[i][1])
        replaceeqs[s] = e
    neweqsdict = dict((newspecies[s], neweqs[s]) for s in range(len(newspecies)))

    if debug:
        for s in newspecies: print("{}: {}".format(s, replaceeqs[s] - neweqsdict[s]).factor())

    return sum((replaceeqs[s] - neweqsdict[s]).factor() != 0 for s in newspecies)


def adair_two_sites():
    """Adair equation for a protein that binds ligand at two identical sites (Ingalls 3.3.1)."""
    print("Adair equation.")
    crn = from_react_file(os.path.join(input_reactions, "adair_two_sites"))
    for p in ['PX2', 'PX1']: crn.qss(p)
    saturation = sympify('(PX1 + 2 * PX2)/2/(P + PX1 + PX2)')
    for s, expr in crn.removed_species: saturation = saturation.subs(s, expr).factor()
    Y = sympify("(X/K1 + X**2/(K1*K2))/(1 + 2*X/K1 + X**2/(K1*K2))").subs(sympify("K1"), sympify("k_1/k1")). \
                            subs(sympify("K2"), sympify("k_2/k2"))
    fail_if_not_equal(sp.factor(saturation - Y), 0)


def allosteric_activation():
    """Allosteric activation (Ingalls 3.7.8)."""
    print("Allosteric activation.")
    crn = from_react_file(os.path.join(input_reactions, "allosteric_activation"))
    crn.qss(cons_law = ('E', ConsLaw('E + ER + ERS', 'Etot')), \
            remove_const = True, merge_reacts = True)
    #crn.remove_all_constants()
    inds = crn.complexes.index(parse_complex('S'))
    fail_if_not_equal(sp.factor(crn.laplacian[inds,inds] - sympify("R*k3*Etot/(R * (k_2 + k3)/k2 + k_1*(k_2 + k3)/(k1*k2) + S*R)")), 0)


def atp_substrate_inhibitor():
    """ATP acts as substrate and inhibitor."""
    print("ATP example.")
    crn = from_react_file(os.path.join(input_reactions, "atp_substrate_inhibitor"))

    # Removing eatpi and eatp using conservation law
    crn.remove(qss = ['eatpi', 'eatp'],
               cons_law = ('e', ConsLaw('e + eatp + eatpi', 'et')))
    indatp = crn.complexes.index(parse_complex('atp'))
    fail_if_not_equal(sp.factor(crn.laplacian[indatp,indatp] - sympify("(k2 * et)/((k_1 + k2) / k1 + atp + atp**2 / (k_3 / k3))")), 0)

    # Saving and reading
    crn.save_reaction_file(os.path.join(input_reactions, "atp_substrate_inhibitor_simplified"))
    crn = from_react_file(os.path.join(input_reactions, "atp_substrate_inhibitor_simplified"))

    # Removing eatpi and eatp without using conservation law
    crn = from_react_file(os.path.join(input_reactions, "atp_substrate_inhibitor"))
    crn.qss('eatpi')
    crn.qss('eatp')
    fail_if_not_equal(sp.factor(crn.laplacian[0,0] - sympify("(k2 * k1)/(k_1 + k2)")), 0)


def competitive_inhibition():
    """Competitive inhibition (Ingalls 3.2.1)."""
    print("Competitive inhibition.")
    crn = from_react_file(os.path.join(input_reactions, "competitive_inhibition"))
    # Using conservation law
    crn.qss(cons_law = ('E', ConsLaw('E + C + C_i', 'Et')))
    fail_if_not_equal(sp.factor(crn.laplacian[0,0] - sympify("k2*Et/(I*((k_1+k2)/k1)/(k_3/k3) + S + ((k_1+k2)/k1))")), 0)

    # Without using conservation law
    crn = from_react_file(os.path.join(input_reactions, "competitive_inhibition"))
    crn.qss('C_i')
    crn.qss('C')
    crn.remove_all_constants()


def double_displ():
    """Double displacement or ping pong mechanism (Ingalls 3.7.6)."""
    print("Double displacement.")
    crn = from_react_file(os.path.join(input_reactions, "double_displacement"))
    origeqs = crn.equations()
    origspecies = crn.species
    crn.qss(cons_law = ('e' , ConsLaw('e + c1 + c2 + ee', 'et')))
    fail_if_not_equal(eqs_match(origeqs, origspecies, crn.removed_species, crn.equations(), crn.species), 0)


def enzyme_kinetics():
    """One-substrate enzyme kinetics."""
    print("One-substrate enzyme kinetics.")

    filename = os.path.join(input_sbml, "enzyme.xml")
    crn = from_sbml(filename)
    rate = sympify("_comp*E*vcat_kcat*veq_kon*S/(vcat_kcat + veq_koff)")
    crn.qss('ES')
    fail_if_not_equal((crn.rates[0] - rate).factor(), 0)
    crn.remove_constant('E')
    fail_if_not_equal((crn.rates[0] - rate).factor(), 0)

    # Michaelis-Menten
    crn = from_sbml(filename)
    enzyme_cons_law = ConsLaw('E + ES', 'Et')
    crn.qss(cons_law = ('E', enzyme_cons_law))
    rate = sympify("_comp*Et*vcat_kcat*veq_kon*S/(vcat_kcat + veq_koff + veq_kon*S)")
    fail_if_not_equal((crn.rates[0] - rate).factor(), 0)

    # Rapid equilibrium
    filename = os.path.join(input_sbml, "enzyme.xml")
    crn = from_sbml(filename)
    enzyme_cons_law = ConsLaw('E + ES', 'Et')
    crn.rapid_eq('ES', 'S + E', cons_law = ('E', enzyme_cons_law))
    fail_if_not_equal((crn.kinetic_params[0] - sympify("Et*vcat_kcat*_comp/(S + veq_koff/veq_kon)")).factor(), 0)


def enzyme_reversible():
    """One-substrate enzyme reversible kinetics."""
    print("One-substrate reversible enzyme kinetics.")
    crn = from_react_file(os.path.join(input_reactions, "enzyme_reversible"))

    crn.qss(cons_law = ('E', ConsLaw('E + C', 'Et')))
    forward = sympify("Vf*S/Ks/(1 + S/Ks + P/Kp)").subs("Vf", sympify("Et*k2"))\
                                                  .subs("Ks", sympify("(k_1+k2)/k1"))\
                                                  .subs("Kp", sympify("(k_1+k2)/k_2")).factor()
    backward = sympify("Vr*P/Kp/(1 + S/Ks + P/Kp)").subs("Vr", sympify("Et*k_1"))\
                                                   .subs("Ks", sympify("(k_1+k2)/k1"))\
                                                   .subs("Kp", sympify("(k_1+k2)/k_2")).factor()
    fail_if_not_equal(set(crn.rates), set([forward, backward]))


def non_competitive_inhibition():
    """Non-competitive inhibition (Segel, Enzyme kinetics)."""
    print("Non-competitive inhibition.")
    crn = from_react_file(os.path.join(input_reactions, "noncompetitive_inhibition"))
    crn.remove(rapid_eq = [('ei', 'e + i'), ('esi', 'e + s + i'), ('es', 's + e')], \
               cons_law = ('e', ConsLaw('e + ei + es + esi', 'et')))
    fail_if_not_equal(sp.factor(crn.laplacian[0,0] - sympify("et*k2/((1 + k3/k_3 * i)*(s + k_1/k1))")), 0)


def passive_transport():
    """Passive transport (Ingalls, section 3.4.2) """
    print("Passive transport.")
    crn = from_react_file(os.path.join(input_reactions, "passive_transport_simpl"))
    crn.qss(cons_law = ('T', ConsLaw('T + TS', 'Ttot')))
    forward = sympify("a1*S1/K1/(1 + S1/K1 + S2/K2)").subs("a1", sympify("k2*Ttot"))\
                                                     .subs("K1", sympify("(k_1+k2)/k1"))\
                                                     .subs("K2", sympify("(k_1+k2)/k_2")).factor()
    backward = sympify("a2*S2/K2/(1 + S1/K1 + S2/K2)").subs("a2", sympify("k_1*Ttot"))\
                                                      .subs("K1", sympify("(k_1+k2)/k1"))\
                                                      .subs("K2", sympify("(k_1+k2)/k_2")).factor()
    fail_if_not_equal(set(crn.rates), set([forward, backward]))

    # Rapid equilibrium on transport step
    crn = from_react_file(os.path.join(input_reactions, "passive_transport"))
    crn.rapid_eq_with_pool('TS1', 'TS2', pool_name = "TS")

    # Assuming only qqs
    crn = from_react_file(os.path.join(input_reactions, "passive_transport"))
    crn.qss(cons_law = ('T', ConsLaw('T + TS1 + TS2', 'Ttot')))


def ternary_compulsory():
    """Transfer of a radioactive atom (a_s -> p_s), in a ternary compulsory mechanism. Cornish-Bowden, 6.8."""
    print("Example of transfer of a radioactive atom.")
    # Exchange requires a_s to bind to e. Inhibited by high concentrations of a and q, as they also bind to e.
    crn = from_react_file(os.path.join(input_reactions, "ternary_compulsory"))
    with warnings.catch_warnings(record = True) as w:
        # Use 'e + ea + eab + eq' instead of proper conservation 'e + ea + ea_s + eab + ea_sb + eq')
        crn.remove(rapid_eq = [('eab', 'ea + b'), ('ea', 'e + a'), ('eq', 'e + q')], \
                             cons_law = ('e' , ConsLaw('e + ea + eab + eq', 'et')))

        crn.remove(qss = ['ea_sb', 'ea_s'])
        # unlabelled species assumed constant
        for s in ['a', 'p', 'b', 'q']:
            crn.remove_constant(s)
        for ws in w: assert "not constant." in str(ws.message)

    # page 122
    ratea_s = sympify('k1*k2*k3*et*b/((1+k1*a/k_1+k1*k2*a*b/(k_1*k_2)+k_4*q/k4)*(k_1*(k_2+k3)+k2*k3*b))')
    ratep_s = sympify('k_1*k_2*k_3*k_4/k4*et*q/((1+k1*a/k_1+k1*k2*a*b/(k_1*k_2)+k_4*q/k4)*(k_1*(k_2+k3)+k2*k3*b))')
    ratea_s = ratea_s.expand().factor().factor()
    ratep_s = ratep_s.expand().factor().factor()
    inda_s = crn.complexes.index(parse_complex('a_s'))
    indp_s = crn.complexes.index(parse_complex('p_s'))
    diffa_s = crn.laplacian[inda_s, inda_s] - ratea_s
    diffa_s = diffa_s.factor()
    diffp_s = crn.laplacian[indp_s, indp_s] - ratep_s
    diffp_s = diffp_s.factor()
    fail_if_not_equal(diffa_s, 0)
    fail_if_not_equal(diffp_s, 0)

    # Full version
    crn = from_react_file(os.path.join(input_reactions, "ternary_compulsory"))
    #for c in ['a', 'p', 'b', 'q']: crn.remove_constant(c)
    crn.remove(rapid_eq = [('eab', 'ea + b'), ('ea', 'e + a'), ('eq', 'e + q')], \
                         qss = ['ea_sb', 'ea_s'], \
                         cons_law = ('e' , ConsLaw('e + ea + ea_s + eab + ea_sb + eq', 'et')))


def three_subs_one_prod_compulsory_irr():
    """Cornish-Bowden, example section 6.9, irreversible version."""
    print("Three substrate, one product, compusory irreversible mechanism.")
    crn = from_react_file(os.path.join(input_reactions, "three_subs_mech_irr"))
    origeqs = crn.equations()
    origspecies = crn.species
    fail_if_not_equal(crn.is_constant('e + e1 + ea + er + eab + eqr'), True)
    crn.qss(cons_law = ('e' , ConsLaw('e + e1 + ea + er + eab + eqr', 'et')))
    fail_if_not_equal(eqs_match(origeqs, origspecies, crn.removed_species, crn.equations(), crn.species), 0)
    fail_if_not_equal(set(get_monoms(crn.rates[0].as_numer_denom()[1], ["a", "b", "c"])), \
                     set(map(sympify, ["c", "a*b", "a*c", "b*c", "a*b*c"])))


def three_subs_one_prod_compulsory_irr_0():
    """Irreversible three substrate, one product compulsory order mechanism."""
    print("Three substrate, one product, compusory irreversible mechanism, v2.")
    crn = from_react_file(os.path.join(input_reactions, "three_subs_one_prod_compul_irr"))
    origeqs = crn.equations()
    origspecies = crn.species
    crn.qss(cons_law = ('e' , ConsLaw('e + ea + eab + eabc', 'et')))
    fail_if_not_equal(eqs_match(origeqs, origspecies, crn.removed_species, crn.equations(), crn.species), 0)


def three_subs_one_prod_full():
    """Three substrates, one product, with rapid eq. and qss."""
    print("Three substrate, one product, compusory irreversible mechanism, v3.")
    crn = from_react_file(os.path.join(input_reactions, "three_subs_one_prod_full"))
    origeqs = crn.equations()
    origspecies = crn.species
    fail_if_not_equal(crn.is_constant('e + ea + eb + ec + eab + eac + ebc + eabc'), True)
    crn.remove(rapid_eq = [('ea', 'e + a'), ('eb', 'e + b'), ('ec', 'e + c'), \
                          ('eab', 'e + a + b'), ('eac', 'e + a + c'), ('ebc', 'e + b +c')], \
               qss = ['eabc'], \
               cons_law = ('e' , ConsLaw('e + ea + eb + ec + eab + eac + ebc + eabc', 'et')))
    fail_if_not_equal(set(get_monoms(crn.rates[0].as_numer_denom()[1], ["a", "b", "c"])), \
                     set(map(sympify, ["1", "a", "b", "c", "a*b", "a*c", "b*c", "a*b*c"])))


def two_catalytic_sites():
    """Enzyme with two catalytic sites (Ingalls, exercise 3.3.4)."""
    print("Enzyme with two catalytic sites.")
    crn = from_react_file(os.path.join(input_reactions, "two_catalytic_sites"))
    crn.qss('c', cons_law = ('e', ConsLaw('e + c', 'et')))
    fail_if_not_equal((crn.rates[0] - sympify('k2*et*s**2/((k_1 + k2) / k1 + s**2)')).factor(), 0)


def two_subs_one_prod_compulsory_irr():
    """Irreversible two substrate, one product compulsory order mechanism.
    Rates compared to http://www.cogsys.cs.uni-tuebingen.de/software/SBMLsqueezer/doc/KineticLaws2.pdf."""
    print("Two substrates, one product compulsory irreversible mechanism.")
    # Irreversible
    crn = from_react_file(os.path.join(input_reactions, "two_subs_one_prod_compul_irr"))
    crn.qss(cons_law = ('e' , ConsLaw('e + ea + eab', 'et')))
    rateab = sympify("k3*et/(kia*Kmb)/(1+a/kia+Kma*b/(kia*Kmb)+a*b/(Kmb*kia))").subs(sympify("kia"), "k_1/k1") \
                                                                                  .subs(sympify("Kma"), "k3/k1") \
                                                                                  .subs(sympify("Kmb"), "(k_2+k3)/k2")
    indab = crn.complexes.index(parse_complex('a + b'))
    fail_if_not_equal((rateab - crn.laplacian[indab,indab]).factor(), 0)


def two_subs_one_prod_compulsory_rev():
    """Reversible two substrate, one product compulsory order mechanism.
    Rates compared to http://www.cogsys.cs.uni-tuebingen.de/software/SBMLsqueezer/doc/KineticLaws2.pdf."""
    print("Two substrates, one product compulsory reversible mechanism.")
    crn = from_react_file(os.path.join(input_reactions, "two_subs_one_prod_compul_rev"))
    crn.qss(cons_law = ('e' , ConsLaw('e + ea + eab', 'et')))

    constants = dict(kpluscat = "k3",
                     kminuscat = "k_1*k_2/(k_1+k_2)",
                     kia = "k_1/k1",
                     kip = "k3/k_3",
                     Kma = "k3/k1",
                     Kmb = "(k_2+k3)/k2",
                     Kmp = "k_1*(k_2+k3)/(k_3*(k_1+k_2))")

    rateab = sympify("kpluscat*et/(kia*Kmb)/(1+a/kia+Kma*b/(kia*Kmb)+a*b/(Kmb*kia)+Kma*b*p/(kia*Kmb*kip)+p/Kmp)").subs(constants)
    indab = crn.complexes.index(parse_complex('a + b'))
    diffab = (rateab - crn.laplacian[indab,indab]).factor()
    ratep = sympify("kminuscat*et/Kmp/(1+a/kia+Kma*b/(kia*Kmb)+a*b/(Kmb*kia)+Kma*b*p/(kia*Kmb*kip)+p/Kmp)").subs(constants)
    indp = crn.complexes.index(parse_complex('p'))
    diffp = (ratep - crn.laplacian[indp,indp]).factor()
    fail_if_not_equal(diffab, 0)
    fail_if_not_equal(diffp, 0)


def two_subs_one_prod_random_irr():
    """Irreversible two substrates, one product random order mechanism.
    (random order bi-uni mechanism, http://www.ebi.ac.uk/sbo/main/SBO:0000432)."""
    print("Two substrates, one product random irreversible mechanism.")

    crn = from_react_file(os.path.join(input_reactions, "two_subs_one_prod_rand_irr"))

    crn.remove(rapid_eq = [('ea', 'e+a'), ('eb', 'e+b')], \
               qss = ['eab'], \
               cons_law = ('e', ConsLaw('e + ea + eb + eab', 'et')))

    constants = dict(Vf = sympify("et*k5"),
                     Vr = sympify("et*(k_3+k_4)"),
                     kia = sympify("k_1/k1"),
                     kib = sympify("k_2/k2"),
                     Kmb = sympify("k1*k_2*(k5 + k_3 + k_4)/(k1*k3*k_2 + k2*k4*k_1)"),
                     Kmp = sympify("(k5 + k_3 + k_4)/k_5"))
    rateab = sympify("Vf/(kia*Kmb)/(1+a/kia+b/kib+a*b/(kia*Kmb))").subs(constants)

    indab = crn.complexes.index(parse_complex('a + b'))
    diff = crn.laplacian[indab, indab] - rateab
    diff = diff.factor()
    fail_if_not_equal(diff, 0)


def two_subs_one_prod_random_rev():
    """Reversible two substrates, one product random order mechanism."""
    print("Two substrates, one product random reversible mechanism.")

    crn = from_react_file(os.path.join(input_reactions, "two_subs_one_prod_rand_rev"))
    crn.remove(rapid_eq = [('ea', 'e+a'), ('eb', 'e+b')], \
                         qss = ['eab'], \
                         cons_law = ('e', ConsLaw('e + ea + eb + eab', 'et')))

    constants = dict(Vf = sympify("et*k5"),
                     Vr = sympify("et*(k_3+k_4)"),
                     kia = sympify("k_1/k1"),
                     kib = sympify("k_2/k2"),
                     Kmb = sympify("k1*k_2*(k5 + k_3 + k_4)/(k1*k3*k_2 + k2*k4*k_1)"),
                     Kmp = sympify("(k5 + k_3 + k_4)/k_5"))

    rateab = sympify("Vf/(kia*Kmb)/(1+a/kia+b/kib+a*b/(kia*Kmb)+p/Kmp)").subs(constants)
    ratep = sympify("Vr/Kmp/(1+a/kia+b/kib+a*b/(kia*Kmb)+p/Kmp)").subs(constants)

    rateab = rateab.expand().factor().factor()
    ratep = ratep.expand().factor().factor()
    indab = crn.complexes.index(parse_complex('a + b'))
    indp = crn.complexes.index(parse_complex('p'))
    diffab = crn.laplacian[indab, indab] - rateab
    diffab = diffab.factor()
    diffp = crn.laplacian[indp, indp] - ratep
    diffp = diffp.factor()
    fail_if_not_equal(diffab, 0)
    fail_if_not_equal(diffp, 0)


def two_subs_two_prods_compulsory():
    """Irreversible and reversible two substrates, two products compulsory order mechanism."""
    print("Two substrates, two products compulsory mechanism.")
    # Full version, without using conservation law
    crn = from_react_file(os.path.join(input_reactions, "two_subs_two_prods_compul"))
    crn.qss('eab')
    crn.qss('ea')
    crn.remove_all_constants()

    # Full version
    crn = from_react_file(os.path.join(input_reactions, "two_subs_two_prods_compul"))
    crn.qss(cons_law = ('e' , ConsLaw('e + ea + eab', 'et')))
    rate = sympify("k3*et*a*b/(k_1*(k_2+k3)/(k1*k2)+(k_2+k3)/k2*a+k3/k1*b+a*b)")
    fail_if_not_equal((crn.rates[0] - rate).factor(), 0)

    # Cornish-Bowden version, 6.3
    crn = from_react_file(os.path.join(input_reactions, "two_subs_two_prods_compul_ternary"))
    crn.qss(cons_law = ('e' , ConsLaw('e + ea + eq + eab', 'et')))

    constants = dict(Vf = sympify("k3*k4*et/(k3+k4)"),
                     Vr = sympify("k_1*k_2*et/(k_1+k_2)"),
                     kia = sympify("k_1/k1"),
                     kib = sympify("(k_1+k_2)/k2"),
                     kip = sympify("(k3+k4)/k_3"),
                     kiq = sympify("k4/k_4"),
                     Kma = sympify("k3*k4/k1/(k3+k4)"),
                     Kmb = sympify("(k_2+k3)*k4/k2/(k3+k4)"),
                     Kmp = sympify("k_1*(k_2+k3)/(k_1+k_2)/k_3"),
                     Kmq = sympify("k_1*k_2/(k_1+k_2)/k_4"))

    rateab = sympify("Vf/(kia*Kmb)/(1+a/kia+b*Kma/(kia*Kmb)+p*Kmq/(kiq*Kmp)+q/kiq+a*b/(kia*Kmb)+Kmq*a*p/(kia*Kmp*kiq)+\
                                       Kma*b*q/(kia*Kmb*kiq)+p*q/(Kmp*kiq)+a*b*p/(kia*Kmb*kip)+b*p*q/(kib*Kmp*kiq))").subs(constants)
    ratepq = sympify("Vr/(kiq*Kmp)/(1+a/kia+b*Kma/(kia*Kmb)+p*Kmq/(kiq*Kmp)+q/kiq+a*b/(kia*Kmb)+Kmq*a*p/(kia*Kmp*kiq)+\
                                       Kma*b*q/(kia*Kmb*kiq)+p*q/(Kmp*kiq)+a*b*p/(kia*Kmb*kip)+b*p*q/(kib*Kmp*kiq))").subs(constants)

    rateab = rateab.factor()
    ratepq = ratepq.factor()
    indab = crn.complexes.index(parse_complex('a + b'))
    indpq = crn.complexes.index(parse_complex('p + q'))
    diffab = crn.laplacian[indab, indab] - rateab
    diffab = diffab.factor()
    diffpq = crn.laplacian[indpq, indpq] - ratepq
    diffpq = diffpq.factor()
    fail_if_not_equal(diffab, 0)
    fail_if_not_equal(diffpq, 0)


def two_subs_two_prods_random():
    """Reversible two substrates, two products random order mechanism."""
    print("Two substrates, two products random mechanism.")
    # Cornish-Bowden version, section 6.1
    crn = from_react_file(os.path.join(input_reactions, "two_subs_two_prods_rand"))
    crn.remove(rapid_eq = [('ea', 'e+a'), ('ep', 'e+p'), ('eb', 'e+b'), ('eq', 'e+q')], \
               qss = ['eab', 'epq'], \
               cons_law = ('e', ConsLaw('e + ea + eb + ep + eq + eab + epq', 'et')), \
               merge_reacts = True)

    constants = dict(Vf = sympify("et*k9*(k5 + k6)/(k5 + k6 + k9 + k_9)"),
                     Vr = sympify("et*k_9*(k_3 + k_4)/(k9 + k_3 + k_4 + k_9)"),
                     kia = sympify("k_1/k1"),
                     kib = sympify("k_2/k2"),
                     kip = sympify("k7/k_7"),
                     kiq = sympify("k8/k_8"),
                     Kmb = sympify("k1*k_2*(k5*k9 + k5*k_3 + k5*k_4 + k6*k9 + k6*k_3 + k6*k_4 + k_3*k_9 + k_4*k_9)/((k1*k3*k_2 + k2*k4*k_1)*(k5 + k6 + k9 + k_9))"),
                     Kmp = sympify("k7*k_8*(k5*k9 + k5*k_3 + k5*k_4 + k6*k9 + k6*k_3 + k6*k_4 + k_3*k_9 + k_4*k_9)/((k8*k_6*k_7 + k7*k_5*k_8)*(k9 + k_3 + k_4 + k_9))"))

    rateab = sympify("Vf/(kia*Kmb)/(1+a/kia+b/kib+p/kip+q/kiq+a*b/(kia*Kmb)+p*q/(Kmp*kiq))").subs(constants)
    ratepq = sympify("Vr/(Kmp*kiq)/(1+a/kia+b/kib+p/kip+q/kiq+a*b/(kia*Kmb)+p*q/(Kmp*kiq))").subs(constants)

    indab = crn.complexes.index(parse_complex('a + b'))
    indpq = crn.complexes.index(parse_complex('p + q'))
    diffab = crn.laplacian[indab, indab] - rateab
    diffab = diffab.factor()
    diffpq = crn.laplacian[indpq, indpq] - ratepq
    diffpq = diffpq.factor()
    fail_if_not_equal(diffab, 0)
    fail_if_not_equal(diffpq, 0)


    # Simplified version
    crn = from_react_file(os.path.join(input_reactions, "two_subs_two_prods_rand_v2"))
    crn.remove(rapid_eq = [('ea', 'e+a'), ('ep', 'e+p'), ('eb', 'e+b'), ('eq', 'e+q')], \
                         qss = ['eab'], \
                         cons_law = ('e', ConsLaw('e + ea + eb + ep + eq + eab', 'et')))

    constants = dict(Vf = sympify("et*(k5+k6)"),
                     Vr = sympify("et*(k_3+k_4)"),
                     kia = sympify("k_1/k1"),
                     kib = sympify("k_2/k2"),
                     kip = sympify("k7/k_7"),
                     kiq = sympify("k8/k_8"),
                     Kmb = sympify("k1*k_2*(k5 + k6 + k_3 + k_4)/(k1*k3*k_2 + k2*k4*k_1)"),
                     Kmp = sympify("k7*k_8*(k5 + k6 + k_3 + k_4)/(k8*k_6*k_7 + k7*k_5*k_8)"))

    rateab = sympify("Vf/(kia*Kmb)/(1+a/kia+b/kib+p/kip+q/kiq+a*b/(kia*Kmb)+p*q/(Kmp*kiq))").subs(constants)
    ratepq = sympify("Vr/(Kmp*kiq)/(1+a/kia+b/kib+p/kip+q/kiq+a*b/(kia*Kmb)+p*q/(Kmp*kiq))").subs(constants)

    indab = crn.complexes.index(parse_complex('a + b'))
    indpq = crn.complexes.index(parse_complex('p + q'))
    diffab = crn.laplacian[indab, indab] - rateab
    diffab = diffab.factor()
    diffpq = crn.laplacian[indpq, indpq] - ratepq
    diffpq = diffpq.factor()
    fail_if_not_equal(diffab, 0)
    fail_if_not_equal(diffpq, 0)


def two_subs_two_prods_subs_enzyme():
    """Substituted-enzyme (ping-pong) mechanism. From Cornish-Bowden, section 6.1."""
    print("Substituted enzyme mechanism.")
    crn = from_react_file(os.path.join(input_reactions, "two_subs_two_prods_subs_enzyme"))
    crn.qss(cons_law = ('e' , ConsLaw('e + e1 + ea + e1b', 'et')))

    constants = dict(Vf = sympify("k2*k4*et/(k2+k4)"),
                     Vr = sympify("k_1*k_3*et/(k_1+k_3)"),
                     kia = sympify("k_1/k1"),
                     kib = sympify("k_3/k3"),
                     kip = sympify("k2/k_2"),
                     kiq = sympify("k4/k_4"),
                     Kma = sympify("(k_1+k2)*k4/k1/(k2+k4)"),
                     Kmb = sympify("k2/k3*(k_3+k4)/(k2+k4)"),
                     Kmp = sympify("k_3/k_2*(k_1+k2)/(k_1+k_3)"),
                     Kmq = sympify("k_1/k_4*(k_3+k4)/(k_1+k_3)"))

    rateab = sympify("Vf/(kia*Kmb)/(a/kia+b*Kma/(kia*Kmb)+p/kip+q*Kmp/(kip*Kmq)+a*b/(kia*Kmb)+a*p/(kia*kip)+\
                                       Kma*b*q/(kia*Kmb*kiq)+p*q/(Kmq*kip))").subs(constants)

    ratepq = sympify("Vr/(kip*Kmq)/(a/kia+b*Kma/(kia*Kmb)+p/kip+q*Kmp/(kip*Kmq)+a*b/(kia*Kmb)+a*p/(kia*kip)+\
                                       Kma*b*q/(kia*Kmb*kiq)+p*q/(Kmq*kip))").subs(constants)

    rateab = rateab.expand().factor().factor()
    ratepq = ratepq.expand().factor().factor()
    indab = crn.complexes.index(parse_complex('a + b'))
    indpq = crn.complexes.index(parse_complex('p + q'))
    diffab = crn.laplacian[indab, indab] - rateab
    diffab = diffab.factor()
    diffpq = crn.laplacian[indpq, indpq] - ratepq
    diffpq = diffpq.factor()
    fail_if_not_equal(diffab, 0)
    fail_if_not_equal(diffpq, 0)


if __name__ == "__main__":
    adair_two_sites()
    allosteric_activation()
    atp_substrate_inhibitor()
    competitive_inhibition()
    double_displ()
    enzyme_kinetics()
    enzyme_reversible()
    non_competitive_inhibition()
    passive_transport()
    ternary_compulsory()
    three_subs_one_prod_compulsory_irr()
    three_subs_one_prod_compulsory_irr_0()
    three_subs_one_prod_full()
    two_catalytic_sites()
    two_subs_one_prod_compulsory_irr()
    two_subs_one_prod_compulsory_rev()
    two_subs_one_prod_random_irr()
    two_subs_one_prod_random_rev()
    two_subs_two_prods_compulsory()
    two_subs_two_prods_random()
    two_subs_two_prods_subs_enzyme()
