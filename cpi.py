#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json, os, pickle, re, requests, sys, traceback
import matplotlib.pyplot as plt
import numpy as np
from bs4 import BeautifulSoup
# from IPython import get_ipython
from mlFunctions import tic, toc, beep, alarm1, run_from_ipython
from progress.bar import ChargingBar

'''
Project Name:

    cpi -- Chemical-Protein Interaction

    It is recommended to run this in the command line thus:

    $ ipython3 -i -m cpi

'''

'''
Meta-stuff
'''

# For running IPython magic commands (e.g., %matplotlib)
# ipython = get_ipython()

# Use latex in matplotlib
# plt.rc('text', usetex=True)

# Determine if running on Saturn server (Linux OS) or my personal machine (Darwin OS)
OS = sys.platform

# Relative path setup
DIR0 = os.getcwd() +'/'
if OS == 'darwin':
    DIR = '/Users/Herman/Documents/jzhang/cpi/'
elif OS == 'linux':
    DIR = '/home/herman/cpi/'

# Display plots inline and change default figure size
if OS == 'darwin':
    pass
    # ipython.magic("matplotlib")
elif OS == 'linux':
    # matplotlib.use("Agg")
    pass

'''
Notes
'''

'''
################################################################################
##### Create chemical-protein dictionary #######################################
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

def downloadCidSyns(DIR, cpiDic, savePickle=False):
    '''
    It takes less than 10 minutes to download all the synonyms for the STITCH CPI database on residential broadband
    '''

    UserAgent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/11.1 Safari/605.1.15'
    Headers = {'User-Agent': UserAgent}

    cidList = [int(cid[-9:]) for cid in list(cpiDic.keys())]

    chunkSize = 190 # chunkSize is the number of length-9 CIDs (plus comma) that can fit into a URL, minus the approximately 100 other characters for the PUG request to PubChem servers.
    numChunks = np.ceil(len(cidList)/chunkSize)

    # Troubleshooting
    # Comment this, because it's useful and I'll probably need it in the future.
    oldNumChunks = numChunks
    numChunks = int(np.ceil(numChunks * 0.05))
    print('Troubleshooting\nReplaced old numChunks size (%d) with %d\n' % (oldNumChunks, numChunks))

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
    alarm1()

    # Pickle requestResults
    # I think requestResults cannot be pickled because it's a bs4 object. I can turn them into Python string objects or save as JSON, instead of Pickle
    # For example:
    # R = requestResults[0]
    # l = [thing for thing in R.strings]
    # See https://stackoverflow.com/questions/16367104/why-does-this-pickle-reach-maximum-recursion-depth-without-recursion
    while savePickle:
        try:
            pname = DIR + 'requestResults.pickle'
            with open(pname, 'wb') as handle:
                pickle.dump(requestResults, handle, protocol=pickle.HIGHEST_PROTOCOL)
            whileGate = False
        except RecursionError as err:
            recursionLimit = sys.getrecursionlimit()
            msg = 'Maximum recursion depth (%d) exceeded while calling a Python object. Would you like to increase it?\nEnter the NEW RECURSION LIMIT or press \'ENTER\' to cancel. (The system limit is somewhere between 30,000 and 40,000)' % recursionLimit
            Input = input(msg)
            # Turn input into int, if possible
            if Input == '':
                break
            else:
                try:
                    Input = int(Input)
                    print('Setting new recursion limit to %d' % Input)
                    sys.setrecursionlimit(Input)
                except:
                    print('Incorrect input type. Enter an integer value greater than %d to set a new recursion limit, or press \'ENTER\' to cancel' % recursionLimit)

    return requestResults

def makeCidSynsDic(DIR, cpiDic):
    '''
    Creates a dictionary of CID (compound ID) synonyms. Synonyms are acquired by
    using the PUG REST utility from PubChem (NIH). Requires internet connection.
    THe dictionary is in the format
        STITCH_ID : [synonym1, synonym2, ... ]

    INPUT:  DIR, string, the directory to read from and to.
            cpiDic, a dic, the dictionary of CPIs (chemical-protein interactions)
    OUTPUT: cidSynsDic, a dic, the dictionary of CID synonyms.
    '''

    # The pickle check is unnecessary because the object is a beautifulSoup object, which cannot be pickled. It is highly recursive. There is no time wasted between its conception and the cidSynsDic creation (just a few lines of code), so no point in pickling it.
    # Determine if requestResults pickle exists. If found, load and process. Else, download, save, and process.
    if 'requestResults.pickle' in os.listdir(DIR):
        # Make sure it's not empty
        isEmpty = os.stat(DIR + 'requestResults.pickle').st_size == 0
        if isEmpty:
            requestResults = downloadCidSyns(DIR, cpiDic)
        else:
            requestResults = pickle.load( open( DIR + 'requestResults.pickle', 'rb' ) )
        # Make sure lengths match
        isDifferentLengths = len(requestResults) != len(cpiDic)
        # Inform user of mismatching lengths
        if isDifferentLengths:
            isBadInput = True
            while isBadInput:
                text = 'Warning, the lengths of the pickled synonyms (requestResults) and the dictionary of CPIs are of different lengths. This implies data that do not match. Choose your option:\n\
\n\
a). Download synonyms again\n\
b). Use pickled synonyms\n\
c). Exit program\n\
\n\
Enter \'a\', \'b\', or \'c\'\n'
                Input = input(text)
                regex = re.findall('[abc]', Input)
                if len(regex) != 1:
                    isBadInput = True
                else:
                    Choice = regex[0]
                    isBadInput = Choice not in ['a','b','c'];
                ### Separate if-statements
                if isBadInput:
                    print('\nIncorrect input. Choose \'a\', \'b\', or \'c\'')
                else:
                    isBadInput = False
            if Choice == 'a':
                print('You entered \'a\'. Downloading again')
                requestResults = downloadCidSyns(DIR, cpiDic)
            elif Choice == 'b':
                print('You entered \'b\'. Using available synonyms.')
            elif Choice == 'c':
                print('You entered \'c\'. Exiting...')
                sys.exit()
    else:
        requestResults = downloadCidSyns(DIR, cpiDic)


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

'''
################################################################################
##### Troubleshooting ############################################
################################################################################
'''

'''
from importlib import reload

from cpi import DIR, loadPickle
cpiDic = loadPickle(DIR, 'cpiDic')

try:
    reload(cpi)
except NameError:
    import cpi
from cpi import downloadCidSyns
requestResults = downloadCidSyns(DIR, cpiDic)

'''

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
    if 'makeCpiDic' in args:
        try:
            cpiDic = makeCpiDic(DIR)
        except FileNotFoundError:
            print('\nThe file was not found. Check the \'DIR\' variable to see if it is pointing to a valid directory.\n\'DIR\' is pointing to %s.\n' % DIR)
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
    print('Commands executed. Copy and paste the below lines to load results into the interpreter.\n')
    for r in results:
        print(r)

# Boilerplate
if __name__ == '__main__':
  main()
