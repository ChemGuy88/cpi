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
    from importlib import reload
    try:
        reload(cpi)
        reload(cpiaux)
    except NameError:
        import cpi, cpiaux
    from cpi import *
    from cpiaux import *

if False:
    # Test dictionary makers

    # cpiDic = makeCpiDic()
    # cidList = makeCidList(DIR, dataDir, outDir, fname='9606.protein_chemical.links.v5.0.tsv', verbose=1, quickMode=True, quickModeLimit=500)
    cidList = makeCidList(DIR, dataDir, outDir, fname='9606.protein_chemical.links.v5.0.tsv', verbose=1, quickMode=False, quickModeLimit=500)
    chemSynsDic, chemSynsRDic = makeChemSynsDic(DIR, outDir, cidList, verbose=1)
    testDicCompleteness(chemSynsDic, chemSynsRDic, 2)

    protSynsDic, protSynsRDic = makeProtSynsDic(DIR, dataDir='STITCH Data/', fname='9606.protein.aliases.v10.5.txt', verbose=1, quickMode=True, quickModeLimit=500)

    ctopDic, ptocDic, pairsToLinksDic = makeLinksDic(DIR, dataDir='STITCH Data/', fname='9606.actions.v5.0.tsv', verbose=1, quickMode=True, quickModeLimit=500)

if False:
    # Test dictionaries
    chemSynsDic, chemSynsRDic = loadPickle(DIR,'chemSynsDic'), loadPickle(DIR,'chemSynsRDic')
    protSynsDic, protSynsRDic = loadPickle(DIR,'protSynsDic'), loadPickle(DIR,'protSynsRDic')
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

if False:
    cidList = loadPickle(DIR+'9606Pickles/','9606cidList')
    cpiDic = loadPickle(DIR+'9606Pickles/','9606cpiDic')

if False:
    pass

'''
################################################################################
##### Function Prototypes ######################################################
################################################################################
'''

sys.exit()

def makeCpiDic(DIR, dataDir, outDir, fname='9606.protein_chemical.links.v5.0.tsv', verbose=0, quickMode=False, quickModeLimit=.1):
    '''
    The code is a clone of makeCidList, but the if-block inside the for loop is different.
    '''
    # verbose start
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Running prototype for \'makeCpiDic\' function.' % str(timeStamp))
        TicSum = datetime.timedelta(0,0,0)
        Tic = tic()
    fpath = dataDir + fname

    # Get progress bar length
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Getting progress bar length.' % str(timeStamp))
        Tic = tic()
    if verbose > 0:
        count = getNumLines(fpath, verbose = 0)
    if verbose > 0:
        Toc = toc(Tic)
        TicSum += Toc

    # Validate quickModeLimit
    if quickMode:
        quickModeLimit = getQuickModeLimit(quickModeLimit, count)

    # Main block
    with open(fpath, 'r') as f:

        # create generator
        if quickMode:
            lineGenerator = (f.readline() for _ in range(quickModeLimit))
        else:
            lineGenerator = (f.readline() for _ in range(count))
        headers=next(lineGenerator)

        # initialize list
        cidList = []
        cpiDic = {}

        # The loop depends on the file being used.
        # fileOrigin = getFileOrigin(fname) # not implemented
        fileOrigin = 'simpleLinks'

        # Verbose prelude to if block
        if verbose > 0:
            timeStamp = ts()
            print('\n[%s] Splitting lines.' % str(timeStamp))
            Tic = tic()
        if verbose > 0:
            if quickMode:
                bar = ChargingBar('', max = quickModeLimit)
            else:
                bar = ChargingBar('', max = count)
        # If block
        if fileOrigin == 'simpleLinks':
            for line in lineGenerator:
                cid, c2, c3  = line.split()
                try:
                    cpiDic[cid].append(c2)
                except KeyError:
                    cpiDic[cid] = [c2]
                if verbose > 0:
                    bar.next()
        elif fileOrigin == 'detailedLinks':
            for line in lineGenerator:
                cid, c2, c3, c4, c5, c6, c7 = line.split()
                cidList.append(cid)
                if verbose > 0:
                    bar.next()
        elif fileOrigin == 'actions':
            for line in lineGenerator:
                a, b, c3, c4, c5, c6 = line.split()
                if isCid(a):
                    cidList.append(a)
                elif isCid(b):
                    cidList.append(b)
                else:
                    print('This should not be happening!')
                    sys.exit(0)
                if verbose > 0:
                    bar.next()
    if verbose > 0:
        bar.finish()
        Toc = toc(Tic)
        TicSum += Toc

    # Create set from list
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Creating set from list.' % str(timeStamp))
        Tic = tic()
    cidList = list(set(cidList))
    if verbose > 0:
            Toc = toc(Tic)
            TicSum += Toc

    # Verbose finish
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Done running protoype for \'makeCpiDic\' function.\nTotal elapsed time was %s (h:mm:ss)' % (str(timeStamp), str(TicSum)))

    # return statement
    withHeaders = False
    if withHeaders:
        return headers, cpiDic
    else:
        return cpiDic
