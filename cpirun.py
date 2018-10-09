#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Project Name:

    cpi -- Chemical-Protein Interaction

    Workspace for cpi project
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

    # test all contents
    i = 0
    gate = False
    for prot, list in protSynsDic.items():
        for alias, source in list:
            # finish this
            reverseList = protSynsRDic[alias][0]
            for resultProt, resultSource in reverseList:
                result = resultProt == prot
            print('Protein %s has alias %s. Does it map back successfully?\t %s' % (prot, alias, str(result)))
            i += 1
            if i > 20:
                gate = True
                break
        if gate:
            break
        print('')
