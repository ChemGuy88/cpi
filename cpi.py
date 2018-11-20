#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime, json, os, pickle, re, requests, sys, traceback
import matplotlib.pyplot as plt
import numpy as np
from bs4 import BeautifulSoup
from IPython import get_ipython
from mlFunctions import tic, toc, beep, ts, alarm1, run_from_ipython, getNumLines
from progress.bar import ChargingBar
from progress.spinner import Spinner

'''
Project Name
------------

    cpi -- Chemical-Protein Interaction

File name
---------

    cpi.py

File Description
----------------

    Main module for this project.

    It is recommended to run this in the command line thus:

    $ ipython3 -i -m cpi

    https://stackoverflow.com/questions/36414966/fastest-approach-to-read-a-big-ascii-file-into-a-numpy-array
'''

'''
Meta Global Variables
'''

ipython = get_ipython()
# ipython.magic("matplotlib") # Same as %matplotlib, turn on interactive plotting in IPython

# Determine if running on Saturn server (Linux OS) or my personal machine (Darwin OS)
OS = sys.platform

# Relative path setup
DIR0 = os.getcwd() +'/'
if OS == 'darwin':
    DIR = '/Users/Herman/Documents/jzhang/cpi/'
elif OS == 'linux':
    DIR = '/home/herman/cpi/'

# Use latex in matplotlib
# plt.rc('text', usetex=True) # False by default

'''
Notes
'''

'''
################################################################################
##### Functions ################################################################
################################################################################
'''

def isCid(string):
    x = 0
    if len(string) != 12:
        x += 1
    if string[:3] != 'CID':
        x += 1
    if x == 0:
        return True
    if x != 0:
        # print(x)
        return False

def loadPickle(DIR, pickleName):
    '''
    Loads a pickled object

    INPUT:  DIR, string, the directory to read from and to.
            pickleName, a string, the name of the pickle file. Assumes suffix is not provided.
    RESULT: Returns 'DIR pickleName', a pickle object
    '''

    pname = '%s%s.pickle' % (DIR, pickleName)
    with open(pname, 'rb') as handle:
        pickleObj = pickle.load(handle)
    return pickleObj

def getQuickModeLimit(quickModeLimit, maximum):
    '''
    Makes sure the quickModeLimit is not less than the maximum number of iterables.

    If the limit is a float, it assumes it's the percentage of iterables to operate on.

    If the limit is an integer, it assumes its the number of iterables to operate on.

    If no limit is given (is NoneType), then the smaller of 500 iterations or 1% of the maximum will be returned.
    '''
    if isinstance(quickModeLimit, float):
        quickModeLimit = int(np.floor(quickModeLimit * maximum))
        if quickModeLimit >= maximum:
            quickModeLimit = int(np.floor(0.1 * maximum))
            quickModeLimit = min([500,quickModeLimit])
    elif isinstance(quickModeLimit, int):
        if quickModeLimit >= maximum:
            quickModeLimit = int(np.floor(0.1 * maximum))
            quickModeLimit = min([500,quickModeLimit])
    elif isinstance(quickModeLimit, type(None)):
        quickModeLimit = getQuickModeLimit(1, maximum)
    else:
        # warning
        pass

    return quickModeLimit

def loadFile(DIR, dataDir, fname, withHeaders=False, verbose=0, quickMode=False, quickModeLimit = 10):
    '''
    Reads tab-delimited file and returns it as a list of lists. Scroll to bottom of Docstring for loading times.

    INPUT:  DIR, string, the directory where the STITCH datasets folder are located.
            verbose, interger. 0 -> no verbose output, any value greater than 0 will print progress feedback, and including a progress bar (if not in quickMode)
            quickMode, boolean. If true, will run a shorter number of iterations as determined by \'quickModeLimit\'.
            quickModeLimit, interger. The number of lines from the file to read. Default 10.
    OUTPUT: array, a Numpy array

    Loading times:
    ============================================================================
    9606.actions.v5.0.tsv
    ---------------------
        Takes 16:20 (mm:ss) to load all 22 million human chemical-protein interactions.
    -------------------------
    protein.aliases.v10.5.txt
    -------------------------
        Takes 48:21 (mm:ss) to load all 48 million protein aliases.
    ============================================================================

    >>> import psutil # or psutils?
    >>> mem = psutil.virtual_memory()
    >>> THRESHOLD = 100 * 1024 * 1024  # 100MB
    >>> if mem.available <= THRESHOLD:
    ...     print("warning")
    '''
    # verbose start
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Running \'loadFile\' function.' % str(timeStamp))
        TicSum = datetime.timedelta(0,0,0)
        Tic = tic()
    fpath = DIR + dataDir + fname

    # Get progress bar length
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Getting progres bar length.' % str(timeStamp))
    if not quickMode and verbose:
        with open(fpath) as f:
            lines0 = (f.readline().splitlines()[0] for line in f)
            count = sum(1 for line in lines0)
    if verbose > 0:
        Toc = toc(Tic)
        TicSum += Toc

    # Main block
    with open(fpath) as f:

        # Read file
        if verbose > 0:
            timeStamp = ts()
            print('\n[%s] Reading file.' % str(timeStamp))
        if quickMode:
            lines = (f.readline().splitlines()[0] for line in range(quickModeLimit))
        else:
            lines = (f.readline().splitlines()[0] for line in f)
        headers = next(lines)
        if verbose > 0:
            Toc = toc(Tic)
            TicSum += Toc

        # Split lines
        if verbose > 0:
            timeStamp = ts()
            print('\n[%s] Splitting lines' % str(timeStamp))
            Tic = tic()
        if not quickMode and verbose:
            bar = ChargingBar('', max = count)
        lists = []
        for line in lines:
            row = line.split('\t')
            lists.append(row)
            if not quickMode and verbose:
                bar.next()
        if not quickMode and verbose:
            bar.finish()
        if verbose > 0:
            Toc = toc(Tic)
            TicSum += Toc
            timeStamp = ts()
            print('\n[%s] Done running \'loadFile\' function.\nTotal elapsed time was %s (h:mm:ss)' % (str(timeStamp), str(TicSum)))
        if not quickMode and verbose:
            beep()

    # return statement
    if withHeaders:
        return [[headers]]+lists
    else:
        return lists

def makeCpiDic():
    '''
    See cpirun.py
    '''
    pass

def makeCidList(DIR, dataDir, outDir, fname='9606.protein_chemical.links.v5.0.tsv', verbose=0, quickMode=False, quickModeLimit=.1):
    '''
    Docstring
    '''
    # verbose start
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Running \'makeCidList\' function.' % str(timeStamp))
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
                cidList.append(cid)
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
        print('\n[%s] Done running \'makeCidList\' function.\nTotal elapsed time was %s (h:mm:ss)' % (str(timeStamp), str(TicSum)))

    # return statement
    withHeaders = False
    if withHeaders:
        return headers, cidList
    else:
        return cidList

def downloadCidSyns(cidList, alarm=alarm1, verbose=0):
    '''
    Download all synonyms for a given list of PubChem Compound ID (CID) numbers using PubChem's PUG REST utility. It takes less than 10 minutes to download all the synonyms for the STITCH CPI database on residential broadband.

    INPUT:  cidList, a list, the list of CIDs.
            alarm, an object (default: alarm = alarm1). This object will be called at the end of the function, when all the synonyms are finished downloading. By default \'alarm1\' is called, which is a series of high and low pitch tones executed in Bash using the operating system.
    OUTPUT: requestResults, a BeautifulSoup4 object.
    '''
    # verbose start
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Running \'downloadCidSyns\' function.' % str(timeStamp))
        TicSum = datetime.timedelta(0,0,0)
        Tic = tic()

    UserAgent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/11.1 Safari/605.1.15'
    Headers = {'User-Agent': UserAgent}

    cidList = [int(cid[-8:]) for cid in cidList]

    chunkSize = 190 # chunkSize is the number of length-9 CIDs (plus comma) that can fit into a URL, minus the approximately 100 other characters for the PUG request to PubChem servers.
    numChunks = np.ceil(len(cidList)/chunkSize)

    requestResults = []

    # Loop miscellanea
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Downloading CID information.' % str(timeStamp))
        bar = ChargingBar('', max = numChunks)
    i = 0 # iterations, number of server requests
    j1 = 0 # index dummy variable
    Tic = tic()

    # Loop
    for num in np.arange(numChunks):
        j0 = j1 # starting position index
        j1 = int(min( len(cidList), (num+1) * chunkSize)) # ending position index
        tempList = cidList[j0:j1]
        cidChunk = ','.join([str(element) for element in tempList])
        URL = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+cidChunk+'/synonyms/XML'
        R = requests.get(URL, headers = Headers)
        soup = BeautifulSoup(R.text, 'lxml')
        Found = soup.findAll('information')
        requestResults.extend(Found) # append results for processing later

        # PubChem requires that no more than 5 requests be made per second
        i += 1
        Toc = toc(Tic, mute=True)
        if Toc.total_seconds() <= 1:
            if i == 5:
                time.sleep(1)
                i = 0
                Tic = tic()
            elif i < 5:
                pass
            elif i > 5:
                print('This isn\'t supposed to happen!')
        elif Toc.total_seconds() > 1:
            i = 1
            Tic = tic()

        # Progress bar
        if verbose > 0:
            bar.next()
    if verbose > 0:
        bar.finish()
    if verbose > 0:
        Toc = toc(Tic)
        TicSum += Toc
    if verbose > 1:
        alarm()

    # Verbose finish
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Done running \'makeCidList\' function.\nTotal elapsed time was %s (h:mm:ss)' % (str(timeStamp), str(TicSum)))

    return requestResults

def makeChemSynsDic(DIR, outDir, cidList, verbose=0):
    '''
    Creates a dictionary of CID (compound ID) synonyms. Synonyms are acquired by
    using the PUG REST utility from PubChem (NIH). Requires internet connection.
    The dictionary is in the format
        STITCH_ID : [synonym1, synonym2, ... ]

    INPUT:  DIR, string, the directory to read from and to.
            cidList, a list, the list of CIDs.
    OUTPUT: chemSynsDic, a dic, the dictionary of CID synonyms (i.e., CID -> chemical name).
            chemSynsRDic, a dic, the dictionary of chemical name synonyms (i.e., chemical name -> CID).
    '''
    # verbose start

    # Download synonyms
    requestResults = downloadCidSyns(cidList, verbose=verbose-1)

    # Process html results to list of strings
    chemSynsDic = {}
    for result in requestResults:
        cid = result.cid.text
        syns = [syn.text for syn in result.findAll('synonym')]
        chemSynsDic[cid] = syns

    # create set of chemical names
    namesList = []
    for listOfNames in chemSynsDic.values():
        for name in listOfNames:
            namesList.append(name)
    namesList = set(namesList)

    # prime reverse dictionary
    chemSynsRDic = {}
    for namesList in chemSynsDic.values():
        for name in namesList:
            chemSynsRDic[name] = []

    # populate dictionary
    for cid, listOfNames in chemSynsDic.items():
        for name in listOfNames:
            chemSynsRDic[name].append(cid)

    # Pickle dictionaries
    pname = outDir + 'chemSynsDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(chemSynsDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pname = outDir + 'chemSynsRDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(chemSynsRDic, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return chemSynsDic, chemSynsRDic

def makeProtSynsDic(DIR, dataDir='STITCH Data', fname='/9606.protein.aliases.v10.5.txt', verbose=0, quickMode=False, quickModeLimit=250000):
    '''
    Docstring
    '''

    # Verbose start
    if verbose > 0:
        TicSum = datetime.timedelta(0,0,0)
        timeStamp = ts()
        print('\n[%s] Running \'makeProtSynsDic\' function.' % str(timeStamp))

    protsyns = loadFile(DIR=DIR, dataDir=dataDir, fname=fname, withHeaders=False, verbose=verbose, quickMode=quickMode, quickModeLimit = quickModeLimit)
    # v10.5.txt has 48,366,210 lines that are read in 1h 15m.

    # create set of protein names and aliases
    if verbose > 0:
        timeStamp = ts()
        text = '\n[%s] Creating sets of protein names and aliases.' % str(timeStamp)
        print(text)
        count = len(protsyns)
        bar = ChargingBar('', max = count)
        Tic = tic()
    prots = []
    aliases = []
    for row in protsyns:
        name, alias = row[0], row[1]
        prots.append(name)
        aliases.append(alias)
        if verbose > 0:
            bar.next()
    prots = set(prots)
    aliases = set(aliases)
    if verbose > 0:
        bar.finish()
        Toc = toc(Tic)
        TicSum += Toc

    # Prime dictionary of protein name aliases
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Priming protein synonyms dictionary.' % str(timeStamp))
        count = len(prots)
        bar = ChargingBar('', max = count)
        Tic = tic()
    protSynsDic = {}
    for prot in prots:
        protSynsDic[prot] = []
        if verbose > 0:
            bar.next()
    if verbose > 0:
        print('')
        Toc = toc(Tic)
        TicSum += Toc
        bar.finish()

    # Priming reverse look-up dictionary
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Priming protein synonym reverse look-up dictionary.' % str(timeStamp))
        count = len(aliases)
        bar = ChargingBar('', max = count)
        Tic = tic()
    protSynsRDic = {}
    for alias in aliases:
        protSynsRDic[alias] = []
        bar.next()
    if verbose > 0:
        print('')
        Toc = toc(Tic)
        TicSum += Toc
        bar.finish()

    # Populate dictionaries
    # this takes the most time. For all proteins, it takes 1 hour per percentage point.
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Populating protein synonym and alias dictionaries.' % str(timeStamp))
        count = len(protsyns)
        bar = ChargingBar('', max = count)
        Tic = tic()
    for line in protsyns:
        name, alias, source = line[0], line[1], line[2] # line format is : name alias source
        protSynsDic[name].append((alias, source))
        protSynsRDic[alias].append((name, source))
        if verbose > 0:
            bar.next()
    if verbose > 0:
        bar.finish()
        Toc = toc(Tic)
        TicSum += Toc

    # Pickle dictionaries
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Pickling dictionaries.' % str(timeStamp))
        Tic = tic()
    pname = DIR + 'protSynsDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(protSynsDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pname = DIR + 'protSynsRDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(protSynsRDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
    Toc = toc(Tic)
    TicSum += Toc

    # Verbose exit
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Done running \'makeProtSynsDic\' function.\nTotal elapsed time was %s (h:mm:ss)' % (str(timeStamp), str(TicSum)))
    if not quickMode and verbose:
        beep()

    return protSynsDic, protSynsRDic

def makeLinksDic(DIR, dataDir='STITCH Data/', fname='9606.actions.v5.0.tsv', verbose=0, quickMode=False, quickModeLimit=1000):
    '''
    'links' means interactions.

    Inputs (5):
      1. fname,
      2. DIR,
      3. verbose,
      4. quickMode,
      5. quickModeLimit,

    Output (3):
      Saves three dictionary objects as pickles.
      1. ptocDic,
      2. ctopDic,
      3. pairsToLinksDic,
    '''
    # create dictionary of interactions
    TicSum = datetime.timedelta(0,0,0)

    timeStamp = ts()
    print('\n[%s] Running prototype for \'makeLinksDic\' function.' % str(timeStamp))

    # Load file
    listOfLinks = loadFile(DIR, dataDir, fname, withHeaders=False, verbose=verbose, quickMode=quickMode, quickModeLimit = quickModeLimit)

    # Create set of proteins, CIDs, and links
    if verbose > 0:
        timeStamp = ts()
        text = '\n[%s] Creating sets of protein names, CIDs, CID-protein pairs, and interaction types (links).' % str(timeStamp)
        print(text)
        count = len(listOfLinks)
        bar = ChargingBar('', max = count)
        Tic = tic()
    setProts = []
    setCids = []
    setPairs = []
    setLinks = []
    for line in listOfLinks:
        A, B, link, action, bool, score = line[0], line[1], line[2], line[3], line[4], line[5]
        if isCid(A):
            setCids.append(A)
        else:
            setProts.append(A)
        setPairs.append((A,B))
        setLinks.append(link)
        if verbose > 0:
            bar.next()
    setProts = set(setProts)
    setCids = set(setCids)
    setPairs = set(setPairs)
    setLinks = set(setLinks)
    if verbose > 0:
        bar.finish()
        Toc = toc(Tic)
        TicSum += Toc

    # Prime protein-to-CID dictionary
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Priming protein-to-CID dictionary.' % str(timeStamp))
        count = len(setProts)
        bar = ChargingBar('', max = count)
        Tic = tic()
    ptocDic = {}
    for prot in setProts:
        ptocDic[prot] = []
        if verbose > 0:
            bar.next()
    if verbose > 0:
        print('')
        Toc = toc(Tic)
        TicSum += Toc
        bar.finish()

    # Prime CID-to-proteins dictionary
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Priming CID-to-proteins dictionary.' % str(timeStamp))
        count = len(setCids)
        bar = ChargingBar('', max = count)
        Tic = tic()
    ctopDic = {}
    for cid in setCids:
        ctopDic[cid] = []
        bar.next()
    if verbose > 0:
        print('')
        Toc = toc(Tic)
        TicSum += Toc
        bar.finish()

    # Prime (CID-protein pair)-to-(link type) dictionary
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Priming (CID-protein pair)-to-(link type) dictionary.' % str(timeStamp))
        count = len(setPairs)
        bar = ChargingBar('', max = count)
        Tic = tic()
    pairsToLinksDic = {}
    for pair in setPairs:
        pairsToLinksDic[pair] = []
        bar.next()
    if verbose > 0:
        print('')
        Toc = toc(Tic)
        TicSum += Toc
        bar.finish()

    # Populate dictionaries
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Populating protein-to-CIDs and CID-to-proteins dictionaries.' % str(timeStamp))
        count = len(listOfLinks)
        bar = ChargingBar('', max = count)
        Tic = tic()
    for line in listOfLinks:
        A, B, link, action, bool, score = line[0], line[1], line[2], line[3], line[4], line[5]
        if isCid(A):
            ctopDic[A].append(B)
        else:
            ptocDic[A].append(B)
        pairsToLinksDic[(A,B)].append((link, action, bool, score))
        if verbose > 0:
            bar.next()
    if verbose > 0:
        bar.finish()
        Toc = toc(Tic)
        TicSum += Toc

    # Pickle dictionaries
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Pickling dictionaries.' % str(timeStamp))
        Tic = tic()
    pname = DIR + 'ptocDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(ptocDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pname = DIR + 'ctopDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(ctopDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
    pname = DIR + 'pairsToLinksDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(pairsToLinksDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
    Toc = toc(Tic)
    TicSum += Toc

    # Verbose exit
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Done running \'makeProtSynsDic\' function.\nTotal elapsed time was %s (h:mm:ss)' % (str(timeStamp), str(TicSum)))
    if not quickMode and verbose:
        beep()

    # Return statement
    return ctopDic, ptocDic, pairsToLinksDic


'''
################################################################################
##### Workspace ################################################################
################################################################################
'''

# see cpirun.py for workspace

'''
################################################################################
##### Define command line arguments ############################################
################################################################################
'''

def main():
    # Can add DIR option to overwrite default DIR value
    if len(sys.argv) > 1:
        args = sys.argv[1:]
    else:
        msg = '\ncpi.py requires arguments to run. Try the following:\n\n\
        \'makeCpiDic\'          -- create CPI dictionary\n\
        \'loadCpiDic\'          -- load CPI dictionary\n\
        \'makeChemSynsDic\'      -- create CID synonyms dictionary\n\
        \'loadCidSynsDic\'      -- load the CID synonyms dictionary\n\
        \n'
        print(msg)

        # Temporary, until bug resolved
        if run_from_ipython():
            sys.exit(0)
        else:
            url = 'https://stackoverflow.com/questions/48571212/why-is-sys-exit-causing-a-traceback'
            print('Traceback is printed on exit when using \'python -i\'. This is a Python bug that has already been reported. Note this does not happen with \'IPython -i\'. Refer to %s.\n' % url)
            sys.exit(0)

    # Work sequence
    results = []
    # quickMode is good for testing code by shortening loops. It's also good for troubleshooting. That's why I replaced Troubleshooting with quickMode, it's a more accurate description of the variable.
    if '-quickMode' in args:
        quickMode = True
        # Add quickMode argument to each of the make() functions.
    if 'makeCpiDic' in args:
        try:
            cpiDic = makeCpiDic(DIR)
        except FileNotFoundError:
            print('\nmakeCpiDic error.\n\
    The file was not found. Check the \'DIR\' variable to see if it is pointing to a valid directory, or if the STITCH data is present.\n\'DIR\' is pointing to %s.\n' % DIR)
        else:
            results.append('# cpiDic was created.')
    if 'loadCpiDic' in args:
        try:
            cpiDic = loadPickle(DIR, 'cpiDic') # It'd be great if I could load this directly to interpreter.
        except FileNotFoundError:
            print('\nloadCpiDic error.\n\
    The file was not found. Check the \'DIR\' variable to see if it is pointing to a valid directory or if the file exists.\n\'DIR\' is pointing to %s.\n' % DIR)
        except:
            print("Unexpected error:", sys.exc_info()[0])
        else:
            results.append('cpiDic = loadPickle(DIR, \'cpiDic\')')
    if 'makeChemSynsDic' in args:
        try:
            cpiDic = loadPickle(DIR, 'cpiDic')
        except FileNotFoundError:
            print('\nmakeChemSynsDic error.\n\
    The file was not found. Try running \'makeCpiDic\' first, or check the \'DIR\' variable to see if it is pointing to a valid directory.\n\'DIR\' is pointing to %s.\n' % DIR)
        except:
            print("Unexpected error:", sys.exc_info()[0])
        else:
            chemSynsDic = makeChemSynsDic(DIR, cpiDic)
            results.append('# chemSynsDic was created.')
    if 'loadCidSynsDic' in args:
        try:
            chemSynsDic = loadPickle(DIR, 'chemSynsDic')
        except FileNotFoundError:
            print('\nloadCidSynsDic error.\n\
    The file was not found. Check the \'DIR\' variable to see if it is pointing to a valid directory or if the file exists.\n\'DIR\' is pointing to %s.\n' % DIR)
        else:
            results.append('chemSynsDic = loadPickle(DIR, \'chemSynsDic\')')
    # Convenient list to copy and paste into IPython
    if len(results) > 0:
        print('Commands executed. Copy and paste the below lines to load results into the interpreter.\n')
        for r in results:
            print(r)

# Boilerplate
if __name__ == '__main__':
  main()
