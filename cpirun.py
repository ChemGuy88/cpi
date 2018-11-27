#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Project Name
------------

    cpi -- Chemical-Protein Interaction

File name
---------

    cpirun.py


File Description
----------------

    Workspace for cpi project. Do not run module. Copy and paste code to interpreter.
'''

if False:
    # Global variables
    DIR = '/Users/Herman/Documents/jzhang/cpi/'
    dataDir = '/Users/Herman/Documents/jzhang/cpi/STITCH Data/'
    outDir = '/Users/Herman/Documents/jzhang/cpi/'

    # Modules
    from importlib import reload, import_module
    try:
        reload(cpi)
    except NameError:
        import cpi
    try:
        reload(cpiaux)
    except NameError:
        import cpiaux
    from cpi import *
    from cpiaux import *

    #
    # Current workspace
    #
    # cidList = loadPickle(DIR+'9606Pickles/','9606cidList')

if False:
    #
    # Test maker functions
    #
    cpiDic, cpiRDic = makeCpiDic(DIR, dataDir, outDir, fname='9606.protein_chemical.links.v5.0.tsv', verbose=1, quickMode=False, quickModeLimit=.1)

    cidList = makeCidList(DIR, dataDir, outDir, fname='9606.protein_chemical.links.v5.0.tsv', verbose=1, quickMode=False, quickModeLimit=500)

    chemSynsDic, chemSynsRDic = makeChemSynsDic(DIR, outDir, cidList, verbose=1)

    protSynsDic, protSynsRDic = makeProtSynsDic(DIR, dataDir='STITCH Data/', fname='9606.protein.aliases.v10.5.txt', verbose=1, quickMode=True, quickModeLimit=500)

    ctopDic, ptocDic, pairsToLinksDic = makeLinksDic(DIR, dataDir='STITCH Data/', fname='9606.actions.v5.0.tsv', verbose=1, quickMode=True, quickModeLimit=500)

if False:
    # Test dictionaries
    chemSynsDic, chemSynsRDic = loadPickle(DIR,'chemSynsDic'), loadPickle(DIR,'chemSynsRDic')
    protSynsDic, protSynsRDic = loadPickle(DIR,'protSynsDic'), loadPickle(DIR,'protSynsRDic')
    cpiDic = loadPickle(DIR+'9606Pickles/','9606cpiDic') # need cpiRDic
    ctopDic, ptocDic = loadPickle(DIR,'ctopDic'), loadPickle(DIR,'ptocDic')
    testDicCompleteness(chemSynsDic, chemSynsRDic, 0)
    testDicCompleteness(protSynsDic, protSynsRDic, 0)
    testDicCompleteness(ctopDic, ptocDic, 0)

if False:
    # Pickle Dictionaries
    # pname = outDir + '9606Pickles/' + '9606cidList.pickle'
    # with open(pname, 'wb') as handle:
    #     pickle.dump(cidList, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # pname = outDir + '9606Pickles/' + '9606cpiDic.pickle'
    # with open(pname, 'wb') as handle:
    #     pickle.dump(linksDic, handle, protocol=pickle.HIGHEST_PROTOCOL)

'''
################################################################################
##### Function Prototypes ######################################################
################################################################################
'''
