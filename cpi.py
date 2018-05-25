#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, pickle, sys
import matplotlib.pyplot as plt
import numpy as np
from IPython import get_ipython
from mlFunctions import tic, toc, beep, beeps, alarm
from progress.bar import ChargingBar

# Determine if running on Saturn server (Linux OS) or my personal machine (Darwin OS)
OS = sys.platform

# For running IPython magic commands (e.g., %matplotlib)
ipython = get_ipython()

# Display plots inline and change default figure size
if OS == 'darwin':
    ipython.magic("matplotlib")
elif OS == 'linux':
    # matplotlib.use("Agg")
    pass


# Use latex in matplotlib
plt.rc('text', usetex=True)

'''
Notes
'''

'''
cpi -- Chemical-Protein Interaction
'''

'''
Meta-stuff
'''
# Relative path setup
DIR0 = os.getcwd() +'/'
if OS == 'darwin':
    DIR = '/Users/Herman/Documents/jzhang/cpi/'
elif OS == 'linux':
    DIR = '/home/herman/cpi/'

'''
Create chemical-protein dictionary
'''
if False:
    print('\nOpening file')
    Tic = tic()
    fname = DIR + 'STITCH Data/protein_chemical.links.v1.0.tsv'
    f = open(fname)
    lines = f.readlines()
    headers = lines[0]
    lines = lines[1:]
    Toc = toc(Tic)

    # Pre-process data
    print('\nPreprocessing file')
    Tic = tic()
    Keys = []
    lolines = []
    for line in lines:
        l = line.split()
        lolines.append(l)
        Keys.append(l[0])
    Set = set(Keys)
    Toc = toc(Tic)

    # Initialize dictionary
    print('\nInitializing dictionary')
    Tic = tic()
    cpi_dic = {}
    for key in Set:
        cpi_dic[key] = []
    Toc = toc(Tic)

    # Create dictionary
    print('\nCreating dictionary')
    Tic = tic()
    for l in lolines:
        cpi_dic[l[0]].append( (l[1], l[2]))
    Toc = toc(Tic)

    # Pickle dictionary
    pname = DIR + 'cpi_dic.pickle'
    with open(pname, 'wb') as handle:
        pickle.dump(cpi_dic, handle, protocol=pickle.HIGHEST_PROTOCOL)

'''
Load pickled cpi dictionary
'''
# Load pickle
pname = DIR + 'cpi_dic.pickle'
with open(pname, 'rb') as handle:
    cpi_dic = pickle.load(handle)

'''
Create synonyms dictionary
'''
import requests
from bs4 import BeautifulSoup

# Format:
# STITCH_ID : [alias1, alias2, ...]

UserAgent = 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_13_4) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/11.1 Safari/605.1.15'
Headers = { 'User-Agent': UserAgent}

cidList = [int(cid[-9:]) for cid in list(cpi_dic.keys())]

chunkSize = 190 # chunkSize is the number of length-9 CIDs (plus comma) that can fit into a URL, minus the approximately 100 other characters for the PUG request to PubChem servers.
numChunks = np.ceil(len(cidList)/chunkSize)
requestResults = []

# ###
numChunks = 5; print('Remove this scaffolding line\n')

# Loop miscellanea
bar = ChargingBar('Downloading CID information', max = numChunks)
i = 0 # iterations, number of server requests
j1 = 0 # index dummy variable
Tic = tic()

'''
i1 = 0
bar = ChargingBar('Downloading CID information', max = numChunks)
for num in np.arange(numChunks):
    i0 = i1
    i1 = min( len(b), (num+1) * chunkSize )
    print(b[i0:i1])
    time.sleep(1)
    bar.next()
bar.finish()
'''

# Loop
for num in np.arange(numChunks):
    j0 = j1 # starting position index
    j1 = min( len(cidList), (num+1) * chunkSize) # ending position index
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
            continue
        elif i > 5:
            print('This isn\'t supposed to happen!')
    elif Toc.total_seconds() > 1:
        i = 1
        Tic = tic()

    # Progress bar
    bar.next()
bar.finish()

# Process html results to string lists
cidSyns = {}
for cid in requestResults:
    name = cid.cid.text
    syns = [syn.text for syn in cid.findAll('synonym')]
    cidSyns[name] = syns

# Pickle cid_syns
pname = DIR + 'cidSyns.pickle'
with open(pname, 'wb') as handle:
    pickle.dump(cidSyns, handle, protocol=pickle.HIGHEST_PROTOCOL)
