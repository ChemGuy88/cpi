#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import datetime, json, os, pickle, re, requests, sys, traceback
import matplotlib.pyplot as plt
import numpy as np
from bs4 import BeautifulSoup
from IPython import get_ipython
from mlFunctions import tic, toc, beep, alarm1, run_from_ipython
from progress.bar import ChargingBar

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
##### Load tab-delimited text file #############################################
################################################################################
'''

def rowFromString(lists):
    '''
    ...
    '''
    pass

def loadFile(fileName, DIR, withHeaders=False, verbose=0, quickMode=False, quickModeLimit=10):
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

    >>> import psutils
    >>> mem = psutils.virtual_memory()
    >>> THRESHOLD = 100 * 1024 * 1024  # 100MB
    >>> if mem.available <= THRESHOLD:
    ...     print("warning")
    '''
    if verbose > 0:
        print('\nRunning \'loadFile\' function.')
        TicSum = datetime.timedelta(0,0,0)
        Tic = tic()
    fname = DIR + fileName

    # Get progress bar length
    if not quickMode and verbose:
        with open(fname) as f:
            lines0 = (f.readline().splitlines()[0] for line in f)
            count = sum(1 for line in lines0)

    # Main block
    with open(fname) as f:

        # Read file
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
            print('\nSplitting lines')
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
            print('\nDone running \'loadFile\' function.\nTotal elapsed time was %s (h:mm:ss)' % str(TicSum))

    if not quickMode and verbose:
        beep()
    if withHeaders:
        return [[headers]]+lists
    else:
        return lists

'''
################################################################################
##### Create chemical-protein interactions dictionary ##########################
################################################################################
'''

def makeCpiDic(DIR):
    '''
    Creates and pickles a dictionary of Chemical-Protein Interactions

    INPUT:  DIR, string, the directory to read from and to.
    OUTPUT: cpiDic, a dictionary
            cpiDic.pickle, the pickled dictionary in directory [DIR]
    '''
    # print('\nOpening file')
    Tic = tic()
    fname = DIR + 'STITCH Data/protein_chemical.links.v1.0.tsv'
    f = open(fname)
    lines = f.readlines()
    headers = lines[0]
    lines = lines[1:]
    # Toc = toc(Tic)

    # Pre-process data
    # print('\nPreprocessing file')
    Tic = tic()
    Keys = []
    lolines = []
    for line in lines:
        l = line.split()
        lolines.append(l)
        Keys.append(l[0])
    Set = set(Keys)
    # Toc = toc(Tic)

    # Initialize dictionary
    # print('\nInitializing dictionary')
    Tic = tic()
    cpiDic = {}
    for key in Set:
        cpiDic[key] = []
    # Toc = toc(Tic)

    # Create dictionary
    # print('\nCreating dictionary')
    Tic = tic()
    for l in lolines:
        cpiDic[l[0]].append( (l[1], l[2]))
    # Toc = toc(Tic)

    # Pickle dictionary
    pname = DIR + 'cpiDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(cpiDic, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return cpiDic

'''
################################################################################
##### Load pickled dictionary ##################################################
################################################################################
'''

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

'''
################################################################################
##### Create synonyms dictionary ###############################################
################################################################################
'''

def downloadCidSyns(cpiDic, alarm=alarm1, quickMode=False):
    '''
    Download all synonyms for a given list of PubChem Compound ID (CID) numbers using PubChem's PUG REST utility. It takes less than 10 minutes to download all the synonyms for the STITCH CPI database on residential broadband.

    INPUT:  cpiDic, a Python dictionary. The list of CIDs is extract from this. This could be changed later to be instead cidList, to allow or more general use.
            alarm, an object (default: alarm = alarm1). This object will be called at the end of the function, when all the synonyms are finished downloading. By default \'alarm1\' is called, which is a series of high and low pitch tones executed in Bash using the operating system.
    OUTPUT: requestResults, a BeautifulSoup4 object.
    '''

    UserAgent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/11.1 Safari/605.1.15'
    Headers = {'User-Agent': UserAgent}

    # Extract list of CIDs from dictionary
    cidList = [int(cid[-9:]) for cid in list(cpiDic.keys())]

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

def makeCidSynsDic(DIR, cpiDic):
    '''
    Creates a dictionary of CID (compound ID) synonyms. Synonyms are acquired by
    using the PUG REST utility from PubChem (NIH). Requires internet connection.
    The dictionary is in the format
        STITCH_ID : [synonym1, synonym2, ... ]

    INPUT:  DIR, string, the directory to read from and to.
            cpiDic, a dic, the dictionary of CPIs (chemical-protein interactions)
    OUTPUT: cidSynsDic, a dic, the dictionary of CID synonyms.
    '''

    # Download synonyms
    requestResults = downloadCidSyns(cpiDic)

    # Process html results to list of strings
    cidSynsDic = {}
    for cid in requestResults:
        name = cid.cid.text
        syns = [syn.text for syn in cid.findAll('synonym')]
        cidSynsDic[name] = syns

    # Save requestResults as JSON
    file = open(DIR + 'cidSynsDic.JSON', 'w')
    json.dump(cidSynsDic, file)
    file.close()

    # Pickle cidSynsDic
    pname = DIR + 'cidSynsDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(cidSynsDic, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return cidSynsDic

def makeProtSynsDic(DIR, verbose=0, quickMode=False):
    '''
        Docstring

    Takes about 47 to 75 minutes on \'protein.aliases.v10.5.txt\'
    '''
    TicSum = 0

    protsynsfn = 'STITCH Data/protein.aliases.v10.5.txt'
    protsyns = loadFile(protsynsfn, DIR, withHeaders=False, verbose=verbose, quickMode=quickMode, quickModeLimit = 6448828)
    # v10.5.txt has 48,366,210 lines that are read in 1h 15m.
    # We can use quickMode to read 644,882 lines in 10m.

    # create set of protein names
    if verbose > 0:
        print('Creating set of protein names')
        Tic = tic()
    prots = []
    for row in protsyns:
        name = row[0]
        if name not in prots:
            prots.append(name)
    if verbose > 0:
        Toc = toc(Tic)
        TicSum += Toc

    # save protein names
    f = open(DIR+'protNames.txt','w')
    for name in prots:
        f.write(name+'\n')
    f.close()
    # v10.5.txt should have 9,507,839 proteins
    # This takes 1 minute to do

    # Create dictionary of protein name aliases
    if verbose > 0:
        print('Creating dictionary of protein name aliases')
        Tic = tic()
    if not quickMode and verbose:
        count = len(protsyns)
        bar = ChargingBar('', max = count)
    protSynsDic = {}
    for prot in prots:
        protSynsDic[prot] = []
    for line in protsyns:
        row = protsyns.pop(0)
        name, alias, source = row[0], row[1], row[2] # row format is : name alias source
        for prot in prots:
            if prot in name:
                protSynsDic[prot].append((alias, source))
        if not quickMode and verbose:
            bar.next()
    if not quickMode and verbose:
        bar.finish()
    if verbose > 0:
        Toc = toc(Tic)
        TicSum += Toc
        print('\nDone\nTotal elapsed time was %s (h:mm:ss)' % str(TicSum))
    if not quickMode and verbose:
        beep()

    # Pickle Dictionary
    pname = DIR + 'protSynsDic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(protsyns, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # return
    return protSynsDic

'''
################################################################################
##### Workspace ################################################################
################################################################################
'''

if False:
    from importlib import reload

    try:
        reload(cpi)
    except NameError:
        import cpi
    from cpi import *

if True:
    verbose = True
    quickMode = True
    TicSum = 0

    print('\nRunning prototype for \'makeProtSynsDic\' function.')

    protsynsfn = 'STITCH Data/protein.aliases.v10.5.txt'
    protsyns = loadFile(protsynsfn, DIR, withHeaders=False, verbose=verbose, quickMode=quickMode, quickModeLimit = 64488)
    # v10.5.txt has 48,366,210 lines that are read in 1h 15m.
    # We can use quickMode to read 644,882 lines in 10s.

    # create set of protein names
    if verbose > 0:
        text = '\nCreating set of protein names'
        count = len(protsyns)
        bar = ChargingBar(text, max = count)
        Tic = tic()
    prots = []
    for row in protsyns:
        name = row[0]
        if name not in prots:
            prots.append(name)
        if verbose > 0:
            bar.next()
    if verbose > 0:
        bar.finish()
        Toc = toc(Tic)
        TicSum += Toc

    # save protein names
    f = open(DIR+'protNames.txt','w')
    for name in prots:
        f.write(name+'\n')
    f.close()
    # v10.5.txt should have 9,507,839 proteins
    # This takes 1 minute to do

    # Create dictionary of protein name aliases
    if verbose > 0:
        print('Creating dictionary of protein name aliases')
        count = len(protsyns)
        bar = ChargingBar('', max = count)
        Tic = tic()
    protSynsDic = {}
    for prot in prots:
        protSynsDic[prot] = []
    for line in protsyns:
        row = protsyns.pop(0)
        name, alias, source = row[0], row[1], row[2] # row format is : name alias source
        for prot in prots:
            if prot in name:
                protSynsDic[prot].append((alias, source))
        if not quickMode and verbose:
            bar.next()
    if not quickMode and verbose:
        bar.finish()
    if verbose > 0:
        Toc = toc(Tic)
        TicSum += Toc
        print('\nDone creating protein synonyms dictionary.\nTotal elapsed time was %s (h:mm:ss)' % str(TicSum))
    if not quickMode and verbose:
        beep()

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
        \'makeCidSynsDic\'      -- create CID synonyms dictionary\n\
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
    if 'makeCidSynsDic' in args:
        try:
            cpiDic = loadPickle(DIR, 'cpiDic')
        except FileNotFoundError:
            print('\nmakeCidSynsDic error.\n\
    The file was not found. Try running \'makeCpiDic\' first, or check the \'DIR\' variable to see if it is pointing to a valid directory.\n\'DIR\' is pointing to %s.\n' % DIR)
        except:
            print("Unexpected error:", sys.exc_info()[0])
        else:
            cidSynsDic = makeCidSynsDic(DIR, cpiDic)
            results.append('# cidSynsDic was created.')
    if 'loadCidSynsDic' in args:
        try:
            cidSynsDic = loadPickle(DIR, 'cidSynsDic')
        except FileNotFoundError:
            print('\nloadCidSynsDic error.\n\
    The file was not found. Check the \'DIR\' variable to see if it is pointing to a valid directory or if the file exists.\n\'DIR\' is pointing to %s.\n' % DIR)
        else:
            results.append('cidSynsDic = loadPickle(DIR, \'cidSynsDic\')')
    # Convenient list to copy and paste into IPython
    if len(results) > 0:
        print('Commands executed. Copy and paste the below lines to load results into the interpreter.\n')
        for r in results:
            print(r)

# Boilerplate
if __name__ == '__main__':
  main()
