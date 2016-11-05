#!/usr/bin/env python

"""Functions for creating SBML models
from SBML or reaction files."""

import keyword as kw
import libsbml

from .crncomplex import sympify
from .parsereaction import parse_reaction_file

__author__ = "Elisa Tonello"
__copyright__ = "Copyright (c) 2016, Elisa Tonello"
__license__ = "BSD"
__version__ = "0.0.1"


def model_from_reacts(reactions, level = 3, version = 1):
    """Create an SBML model from a list of reactions,
    with one compartment.
    Return model, document and list of species."""

    # Create the list of species from the
    # reactants and products
    species = sorted(list(set([s for r in reactions for s in r.reactant] +
                              [s for r in reactions for s in r.product])))

    # Create an empty SBMLDocument
    # level 3, version 1
    document = libsbml.SBMLDocument(level, version)

    # Create the model
    model = document.createModel()

    # Create compartment
    c1 = model.createCompartment()
    c1.setId('c1')
    c1.setConstant(True)

    # Create the species
    for sid in species:
        s = model.createSpecies()
        s.setId(sid)
        s.setCompartment('c1')
        s.setConstant(False)
        s.setInitialAmount(0)
        s.setBoundaryCondition(False)
        s.setHasOnlySubstanceUnits(False)

    # Create the reactions
    for reaction in reactions:
        r = model.createReaction()
        r.setId(reaction.reactionid)
        r.setReversible(False)
        r.setFast(False)

        for s in sorted(reaction.reactant):
            react_species = r.createReactant()
            react_species.setSpecies(s)
            react_species.setStoichiometry(int(reaction.reactant[s]))
            react_species.setConstant(True)

        for s in sorted(reaction.product):
            prod_species = r.createProduct()
            prod_species.setSpecies(s)
            prod_species.setStoichiometry(int(reaction.product[s]))
            prod_species.setConstant(True)

        math_ast = libsbml.parseL3Formula(str(reaction.rate).replace('**', '^'))
        kinetic_law = r.createKineticLaw()
        kinetic_law.setMath(math_ast)

    promote_params(model, document)
    return model, document, list(species)


def model_from_react_file(filename):
    """Create an SBML model from a reaction file,
    with one compartment.
    Return model, document and list of species."""
    return model_from_reacts(parse_reaction_file(filename))


def replace_reacts(model, document, reactions):
    """Delete the current reactions and replace them with reactions.
    Update the species.
    Return document, model and list of species."""

    species = sorted(list(set([s for r in reactions for s in r.reactant] +
                              [s for r in reactions for s in r.product])))
    current_species = [model.getSpecies(s).getName()
                       if model.getSpecies(s).getName()
                       else model.getSpecies(s).getId() for s in range(model.getNumSpecies())]
    current_species_ids = dict(zip(current_species, [model.getSpecies(s).getId()
                                                     if model.getSpecies(s).getId()
                                                     else model.getSpecies(s).getName()
                                                     for s in range(model.getNumSpecies())]))
    # Delete any unused species
    deleteSpecies = []
    for s in range(model.getNumSpecies()):
        if model.getSpecies(s).getName() not in species and \
           model.getSpecies(s).getId() not in species:
               deleteSpecies.append(s)

    for s in sorted(deleteSpecies, reverse = True):
        model.getSpecies(s).removeFromParentAndDelete()

    # Add any new species
    for sid in species:
        if sid not in current_species:
            current_species.append(sid)
            current_species_ids[sid] = sid
            s = model.createSpecies()
            s.setId(sid)

    # Dict of parameters
    current_param_ids = dict(zip([model.getParameter(s).getName() if model.getParameter(s).getName()
                                                              else model.getParameter(s).getId() for s in range(model.getNumParameters())],
                                 [model.getParameter(s).getId() if model.getParameter(s).getId()
                                                              else model.getParameter(s).getName() for s in range(model.getNumParameters())]))
    current_map = dict(current_species_ids, **current_param_ids)

    # Dict of reaction ids
    react_ids = dict(zip([model.getReaction(s).getName() if model.getReaction(s).getName()
                                                         else model.getReaction(s).getId() for s in range(model.getNumReactions())],
                                 [model.getReaction(s).getId() if model.getReaction(s).getId()
                                                         else model.getReaction(s).getName() for s in range(model.getNumReactions())]))

    # Remove current reactions
    for r in range(model.getNumReactions()-1, -1, -1):
        model.getReaction(r).removeFromParentAndDelete()

    # Add new reactions
    for reaction in reactions:
        r = model.createReaction()
        if reaction.reactionid in react_ids:
            r.setId(react_ids[reaction.reactionid])
        else:
            if reaction.reactionid[0].isdigit():
                reaction.reactionid = "r" + reaction.reactionid
            r.setId(reaction.reactionid)
        r.setName(reaction.reactionid)
        r.setReversible(False)
        r.setFast(False)

        for s in sorted(reaction.reactant):
            react_species = r.createReactant()
            react_species.setSpecies(current_species_ids[s])
            react_species.setStoichiometry(int(reaction.reactant[s]))
            react_species.setConstant(True)

        for s in sorted(reaction.product):
            prod_species = r.createProduct()
            prod_species.setSpecies(current_species_ids[s])
            prod_species.setStoichiometry(int(reaction.product[s]))
            prod_species.setConstant(True)

        rate = reaction.rate
        for s in rate.free_symbols:
            if str(s) in current_map:
                rate = rate.subs(s, sympify(current_map[str(s)]))
        math_ast = libsbml.parseL3Formula(str(rate).replace('**', '^'))
        kinetic_law = r.createKineticLaw()
        kinetic_law.setMath(math_ast)

    return model, document, list(species)


def model_from_sbml(sbmlxml):
    """Read a model from an sbml file.
    Return model and document."""
    # Read the sbml file
    reader = libsbml.SBMLReader()
    document = reader.readSBMLFromFile(sbmlxml)
    if document.getNumErrors() > 0:
        document.printErrors()
        raise SystemExit("Error while reading the SBML file.")

    # Create the model
    model = document.getModel()
    if model == None:
        raise SystemExit("Error retrieving model.")

    promote_params(model, document)
    _adjust(model, document)
    convert_functions(model, document)
    return model, document


def convert_functions(model, document):
    """Replace functions with formulas."""
    config = libsbml.ConversionProperties()
    if config != None:
        config.addOption('expandFunctionDefinitions')
    status = document.convert(config)
    if status != libsbml.LIBSBML_OPERATION_SUCCESS:
        print('Error: function conversion failed:')
        document.printErrors()


def promote_params(model, document):
    """Change parameters to global, so that they are preserved
    when reactions change."""
    convProps = libsbml.ConversionProperties()
    convProps.addOption("promoteLocalParameters", True, "Promotes all Local Parameters to Global ones")
    if (document.convert(convProps) != libsbml.LIBSBML_OPERATION_SUCCESS):
        raise SystemExit("Error promoting local parameters to global.")


def _fix_name(name):
    # replace some characters
    symbols = ["-", "+", " ", ".", "/", "`", "~", "!",
               "(", ")", "[", "]", "{", "}", "<", ">",
               "'", ",", "*", "^", ":", "#", "@", "&"]
    for symbol in symbols:
        name = name.replace(symbol, "_")

    # sympy does not like symbols starting with a number
    if name[0].isdigit():
        name = "_" + name

    # avoid keywords
    kwlist = kw.kwlist + ['source', 're', 'factor', 'EmptySet', 'E1', 'CC', 'ff', 'GF', 'Gt', 'LC', 'lcm', 'comp']
    while name in kwlist:
        name = "_" + name
    return name


def _adjust(model, document):
    """Change dash, /, parenthesis and other symbols to underscore in
    compartments, species, reaction and parameter ids and names,
    to avoid problems with sympy."""

    # species
    for s in range(model.getNumSpecies()):
        original = model.getSpecies(s).getName()
        if original:
            fixed = _fix_name(original)
            if original != fixed:
                model.getSpecies(s).setName(fixed)

        original = model.getSpecies(s).getId()
        if original:
            fixed = _fix_name(original)
            model.getSpecies(s).setId(fixed)

        original = model.getSpecies(s).getCompartment()
        if original:
            fixed = _fix_name(original)
            model.getSpecies(s).setCompartment(fixed)

    # reactions
    for s in range(model.getNumReactions()):
        original = model.getReaction(s).getName()
        if original:
            fixed = _fix_name(original)
            model.getReaction(s).setName(fixed)

        original = model.getReaction(s).getId()
        if original:
            fixed = _fix_name(original)
            if original != fixed:
                model.getReaction(s).setId(fixed)

        kl = model.getReaction(s).getKineticLaw().getMath()
        fix = False
        if kl:
            for i in range(model.getReaction(s).getNumReactants()):
                original = model.getReaction(s).getReactant(i).getSpecies()
                if original:
                    fixed = _fix_name(original)
                    if original != fixed:
                        fix = True
                        kl.replaceArgument(original, libsbml.ASTNode(libsbml.parseL3Formula(fixed)))
                        model.getReaction(s).getReactant(i).setSpecies(fixed)

            for i in range(model.getReaction(s).getNumProducts()):
                original = model.getReaction(s).getProduct(i).getSpecies()
                if original:
                    fixed = _fix_name(original)
                    if original != fixed:
                        fix = True
                        kl.replaceArgument(original, libsbml.ASTNode(libsbml.parseL3Formula(fixed)))
                        model.getReaction(s).getProduct(i).setSpecies(fixed)

            for i in range(model.getReaction(s).getNumModifiers()):
                original = model.getReaction(s).getModifier(i).getSpecies()
                if original:
                    fixed = _fix_name(original)
                    if original != fixed:
                        fix = True
                        kl.replaceArgument(original, libsbml.ASTNode(libsbml.parseL3Formula(fixed)))
                        model.getReaction(s).getModifier(i).setSpecies(fixed)

            if model.getReaction(s).getCompartment():
                original = model.getReaction(s).getCompartment()
                fixed = _fix_name(original)
                if original != fixed:
                    fix = True
                    kl.replaceArgument(original, libsbml.ASTNode(libsbml.parseL3Formula(fixed)))
                    model.getReaction(s).setCompartment(fixed)

            if fix: model.getReaction(s).getKineticLaw().setMath(kl)

    # parameters
    for s in range(model.getNumParameters()):
        original = model.getParameter(s).getName()
        if original:
            fixed = _fix_name(original)
            if original != fixed:
                for r in range(model.getNumReactions()):
                    kl = model.getReaction(r).getKineticLaw().getMath()
                    if kl:
                        kl.replaceArgument(original, libsbml.ASTNode(libsbml.parseL3Formula(fixed)))
                        model.getReaction(r).getKineticLaw().setMath(kl)
                model.getParameter(s).setName(fixed)

        original = model.getParameter(s).getId()
        if original:
            fixed = _fix_name(original)
            if original != fixed:
                for r in range(model.getNumReactions()):
                    kl = model.getReaction(r).getKineticLaw().getMath()
                    if kl:
                        kl.replaceArgument(original, libsbml.ASTNode(libsbml.parseL3Formula(fixed)))
                        model.getReaction(r).getKineticLaw().setMath(kl)
                model.getParameter(s).setId(fixed)

    # compartments
    for s in range(model.getNumCompartments()):
        original = model.getCompartment(s).getName()
        if original:
            fixed = _fix_name(original)
            if original != fixed:
                for r in range(model.getNumReactions()):
                    kl = model.getReaction(r).getKineticLaw().getMath()
                    if kl:
                        kl.replaceArgument(original, libsbml.ASTNode(libsbml.parseL3Formula(fixed)))
                        model.getReaction(r).getKineticLaw().setMath(kl)
                model.getCompartment(s).setName(fixed)

        original = model.getCompartment(s).getId()
        if original:
            fixed = _fix_name(original)
            if original != fixed:
                for r in range(model.getNumReactions()):
                    kl = model.getReaction(r).getKineticLaw().getMath()
                    if kl:
                        kl.replaceArgument(original, libsbml.ASTNode(libsbml.parseL3Formula(fixed)))
                        model.getReaction(r).getKineticLaw().setMath(kl)
                model.getCompartment(s).setId(fixed)
