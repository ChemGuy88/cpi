#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime, json, os, pickle, re, requests, sys, traceback
import matplotlib.pyplot as plt
import numpy as np
from bs4 import BeautifulSoup
from IPython import get_ipython
from mlFunctions import tic, toc, beep, ts, alarm1, run_from_ipython
from progress.bar import ChargingBar
from progress.spinner import Spinner

'''
Project Name:

    cpi -- Chemical-Protein Interaction

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

def makeCidList(DIR, dataDir='STITCH Data/', fname='protein_chemical.links.v1.0.tsv', verbose=0, quickMode=False, quickModeLimit=None):
    '''
    '''
    # verbose start
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Running \'makeCidList\' function.' % str(timeStamp))
        TicSum = datetime.timedelta(0,0,0)
        Tic = tic()
    fpath = DIR + dataDir + fname

    # Get progress bar length
    if verbose > 0:
        timeStamp = ts()
        print('\n[%s] Getting progress bar length.' % str(timeStamp))
    if not quickMode and verbose:
        with open(fpath) as f:
            lines0 = (f.readline().splitlines()[0] for line in f)
            count = sum(1 for line in lines0)
    if verbose > 0:
        Toc = toc(Tic)
        TicSum += Toc

    # Main block
    with open(fpath, 'r') as f:

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
        cidList = []
        for line in lines:
            row = line.split('\t')
            A, B, link, action, bool, score = line[0], line[1], line[2], line[3], line[4], line[5]
            if isCid(A):
                cidList.append(A)
            else:
                cidList.append(B)
            if verbose > 0:
                bar.next()
        cidList = list(set(cidList))
        if not quickMode and verbose:
            bar.finish()
        if verbose > 0:
            Toc = toc(Tic)
            TicSum += Toc
            timeStamp = ts()
            print('\n[%s] Done running \'makeCidList\' function.\nTotal elapsed time was %s (h:mm:ss)' % (str(timeStamp), str(TicSum)))
        if not quickMode and verbose:
            beep()

    # return statement
    return cidList

def downloadCidSyns(cidList, alarm=alarm1, verbose=0, quickMode=False, quickModeLimit=None):
    '''
    Download all synonyms for a given list of PubChem Compound ID (CID) numbers using PubChem's PUG REST utility. It takes less than 10 minutes to download all the synonyms for the STITCH CPI database on residential broadband.

    INPUT:  cidList, a list, the list of CIDs.
            alarm, an object (default: alarm = alarm1). This object will be called at the end of the function, when all the synonyms are finished downloading. By default \'alarm1\' is called, which is a series of high and low pitch tones executed in Bash using the operating system.
    OUTPUT: requestResults, a BeautifulSoup4 object.
    '''

    UserAgent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/11.1 Safari/605.1.15'
    Headers = {'User-Agent': UserAgent}

    chunkSize = 190 # chunkSize is the number of length-9 CIDs (plus comma) that can fit into a URL, minus the approximately 100 other characters for the PUG request to PubChem servers.
    numChunks = np.ceil(len(cidList)/chunkSize)

    # quickMode; Shorten the list
    if quickMode:
        oldNumChunks = numChunks
        numChunks = int(np.ceil(numChunks * 0.05))
        print('quickMode\nReplaced old numChunks size (%d) with %d\n' % (oldNumChunks, numChunks))

    requestResults = []

    # Loop miscellanea
    bar = ChargingBar('Downloading CID information', max = numChunks)
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
        bar.next()
    bar.finish()
    alarm()

    return requestResults

def makeChemSynsDic(DIR, cidList, verbose=0, quickMode=None, quickModeLimit=None):
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

    # Download synonyms
    requestResults = downloadCidSyns(cidList, verbose=0, quickMode=None, quickModeLimit=None)

    # Process html results to list of strings
    chemSynsDic = {}
    for result in requestResults:
        cid = result.result.text
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
    for cid, listOfNames in namesList:
        for name in listOfNames:
            chemSynsRDic[name].append(cid)

    # Pickle chemSynsDic
    pname = DIR + 'chemSynsDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(chemSynsDic, handle, protocol=pickle.HIGHEST_PROTOCOL)

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

def checkDics(DIR, dic1name='protSynsDic', dic2name='ptocDic', elements='keys'):
    '''
    Does set operations on two dictionaries.
    ======================================================================
    Results for reference:
    ======================================================================
    For dictionaries generated from species 9606 (Humans), the following results are obtained:

    linkProts has 0 unique proteins.
    synsProts has 3105 unique proteins.
    They share 17352 (84.82%) proteins.
    ======================================================================

    '''
    # Rename variables to generic names
    dic1 = loadPickle(DIR, dic1name)
    dic2 = loadPickle(DIR, dic2name)
    dic1keys = set(dic1.keys())
    dic2keys = set(dic2.keys())
    print('%s has %d proteins.\n%s has %d %s.' % (dic1name, len(dic1keys), dic2name, len(dic2keys), elements))
    dic1extra = dic1keys - dic2keys
    dic2extra = dic2keys - dic1keys
    common = dic1keys & dic2keys
    total = dic1keys | dic2keys
    a, b, t, u = len(dic1extra), len(dic2extra), len(common), len(total)
    print('%s has %d unique proteins.\n%s has %d unique %s.\nThey share %d (%0.2f%%) of %s.' % (dic1name, a, dic2name, b, elements, t, 100*t/u), elements)

def blocks(file, size=65536):
    '''
    From https://stackoverflow.com/a/9631635/5478086
    '''
    while True:
        b = file.read(size)
        if not b:
            break
        yield b

def getNumLines(fname):
    '''
    From https://stackoverflow.com/a/9631635/5478086
    '''
    with open(fname, 'r', encoding="utf-8", errors='ignore') as f:
        genObj = (bl.count('\n') for bl in blocks(f))
        s = sum(genObj)
        print('%s has %d lines.' % (fname, s))

def getNumLines(fname, verbose=1):
    '''
    docstring
    '''
    # Verbose start
    if verbose > 0:
        TicSum = datetime.timedelta(0,0,0)
        timeStamp = ts()
        print('\n[%s] Counting number of lines in %s.' % (str(timeStamp), fname))
        # spinner = Spinner('')
        Tic = tic()
    with open(fname) as f:
        n = 0
        lines = (f.readline() for line in f)
        for line in lines:
            n += 1
            if verbose > 0:
                pass
                # spinner.next()
    if verbose > 0:
        # spinner.finish()
        Toc = toc(Tic)
        TicSum += Toc

    # verbose exit
    if verbose > 0:
        print('\n[%s] There were %s lines in %s.' % (str(timeStamp), n, fname))

    return n

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
