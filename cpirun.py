#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Project Name:

    cpi -- Chemical-Protein Interaction

    Workspace for cpi project. Do not run module. Copy and paste code to interpreter.
'''

if False:
    from importlib import reload
    try:
        reload(cpi)
    except NameError:
        import cpi
    from cpi import *

    cidList = makeCidList(DIR, dataDir='STITCH Data/', fname='9606.protein_chemical.links.v5.0.tsv')
    chemSynsDic, chemSynsRDic = makeChemSynsDic(DIR, cidList, verbose=1, quickMode=True, quickModeLimit=None)

    protSynsDic, protSynsRDic = makeProtSynsDic(DIR, dataDir='STITCH Data/', fname='9606.protein.aliases.v10.5.txt', verbose=1, quickMode=True, quickModeLimit=250000)

    ctopDic, ptocDic, pairsToLinksDic = makeLinksDic(DIR, dataDir='STITCH Data/', fname='9606.actions.v5.0.tsv', verbose=1, quickMode=True, quickModeLimit=1000)

    # test protSynsDic forward lookup
    prot = list(protSynsDic.keys())[0]
    aliases = protSynsDic[prot]
    print('Protein ID %s has the following aliases:' % prot)
    for alias, source in aliases:
        print('\t'+alias)

    # test protSynsDic backward lookup
    resultProt, resultSource = protSynsRDic[alias][0]
    result = resultProt == prot
    print('Does %s map back to %s?' % (alias, prot))
    print(result)

    # test some contents
    i = 0
    gate = False
    for prot, list in protSynsDic.items():
        for alias, source in list:
            # finish this
            reverseList = protSynsRDic[alias]
            for (resultProt, resultSource) in reverseList:
                result = resultProt == prot
            print('Protein %s has alias %s. Does it map back successfully?\t %s' % (prot, alias, str(result)))
            i += 1
            if i > 20:
                gate = True
                break
        if gate:
            break

    # test all forward lookup contents
    forwardMisses = 0
    numProts = 0
    numAliases = 0
    gate = False
    for prot, list in protSynsDic.items():
        numProts += 1
        for alias, source in list:
            numAliases += 1
            reverseList = protSynsRDic[alias]
            numResults = 0
            for resultProt, resultSource in reverseList:
                result = resultProt == prot
                numResults += result
            if numResults == 0:
                forwardMisses += 1
                Input = input('A reverse lookup failed. Show context? y/n')
                if Input.lower() == 'y':
                    print('While mapping alias %s back to protein %s, the result was instead %s' % (alias, prot, resultProt))
                Input = input('Exit test? y/n')
                if Input.lower() == 'y':
                    gate = True
                    break
            if gate:
                break
        if gate:
            break
    print('%d look-ups failed.' % forwardMisses)
    numEntries = len(protSynsDic)
    print('The dictionary had %d entries. %d were tested.' % (numEntries, numProts))

    # test all reverse lookup contents
    reverseMisses = 0
    numAliases = 0
    numProts = 0
    gate = False
    for alias, list in protSynsRDic.items():
        numAliases += 1
        for prot, source in list:
            numProts += 1
            reverseList = protSynsDic[prot]
            numResults = 0
            for resultAlias, resultSource in reverseList:
                result = resultAlias == alias
                numResults += result
            if numResults == 0:
                reverseMisses += 1
                Input = input('A reverse lookup failed. Show context? y/n')
                if Input.lower() == 'y':
                    print('While mapping prot %s back to alias %s, the result was instead %s' % (prot, alias, resultAlias))
                Input = input('Exit test? y/n')
                if Input.lower() == 'y':
                    gate = True
                    break
            if gate:
                break
        if gate:
            break
    print('%d reverse look-ups failed.' % reverseMisses)
    numREntries = len(protSynsRDic)
    print('The reverse dictionary has %d entries. %d were tested.' % (numREntries, numAliases))

    if (forwardMisses + reverseMisses) == 0:
        print('The tested dictionaries are complete in both the forward and backward directions.')
    else:
        print('Warning. The dictionaries failed some forward or reverse lookups.')
