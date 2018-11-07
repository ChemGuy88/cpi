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

    cpiaux.py

File Description
----------------

    Holds auxiliary files used in debugging and developing cpi.py
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

def getName(element):
    '''
    When checking dictionaries the values of each key might be a tuple or a string. This returns a singleton string value in either case.

    If the dictionary is a forward dictionary it returns the entity's name value. If the dictionary is a reverse look-up dictionary it returns the entity's alias name value.

    In summary:

    Dictionary type | format
    =============================
    forward         | { entity : [alias1, alias2, ...]}
    reverse         | { alias1 : [entity]}
    OR
    forward         | { entity : [(alias1, source1), (alias2, source), ...]}
    reverse         | { (alias1, source1) : [entity]}

    In either case getName returns the same result:

    Dictionary type | getName result
    =============================
    forward         | alias1
    reverse         | entity
    OR
    forward         | alias1
    reverse         | entity
    '''
    if isinstance(element, str):
        return element
    elif isinstance(element, (tuple, list)):
        return element[0]

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
    print('%s has %d %s.\n%s has %d %s.' % (dic1name, len(dic1keys), elements, dic2name, len(dic2keys), elements))
    dic1extra = dic1keys - dic2keys
    dic2extra = dic2keys - dic1keys
    common = dic1keys & dic2keys
    total = dic1keys | dic2keys
    a, b, t, u = len(dic1extra), len(dic2extra), len(common), len(total)
    print('%s has %d unique %s.\n%s has %d unique %s.\nThey share %d (%0.2f%%) of %s.' % (dic1name, a, elements, dic2name, b, elements, t, 100*t/u), elements)

def testDicCompleteness(forwardDic, reverseDic, quickMode=0, quickModeLimit=20):
    '''
    See if all entity aliases in forwardDic map back to entity using the reverseDic. Pictorially:
    forwardDic[entity] --> alias
    reverseDic[alias] --> entity
    '''
    if quickMode == 0:
        testDicCompleteness_full(forwardDic, reverseDic)
    elif quickMode == 1:
        testDicCompleteness_med(forwardDic, reverseDic, quickModeLimit)
    elif quickMode == 2:
        testDicCompleteness_short(forwardDic, reverseDic)

def testDicCompleteness_short(forwardDic, reverseDic):
    # test forwardDic forward lookup
    entity = list(forwardDic.keys())[0]
    aliases = forwardDic[entity]
    print('Entity %s has the following aliases:' % entity)
    for element in aliases:
        alias = getName(element)
        print('\t'+alias)

    # test forwardDic backward lookup
    # resultEntity, resultSource = reverseDic[alias][0]
    resultEntity = getName(reverseDic[alias][0])
    result = resultEntity == entity
    print('Does %s map back to %s?' % (alias, entity))
    print(result)

def testDicCompleteness_med(forwardDic, reverseDic, quickModeLimit=20):
    # test some contents
    i = 0
    gate = False
    for entity, list in forwardDic.items():
        for element in list:
            alias = getName(element)
            # finish this
            reverseList = reverseDic[alias]
            for element in reverseList:
                resultEntity = getName(element)
                result = resultEntity == entity
            print('Entity %s has alias %s. Does it map back successfully?\t %s' % (entity, alias, str(result)))
            i += 1
            if i > quickModeLimit:
                gate = True
                break
        if gate:
            break

def testDicCompleteness_full(forwardDic, reverseDic):
    '''
    Docstring
    '''
    # forwardDic
    numEntities = len(forwardDic)

    # Progress bar
    if True:
        timeStamp = ts()
        print('\n[%s] Testing all forward lookup dictionary contents' % str(timeStamp))
        count = numEntities
        bar = ChargingBar('', max = count)
        Tic = tic()

    # Test all forward lookup dictionary contents
    forwardMisses = 0
    entityCount = 0
    gate = False
    for entity, list in forwardDic.items():
        entityCount += 1
        for element in list:
            alias = getName(element)
            reverseList = reverseDic[alias]
            numResults = 0
            for element in reverseList:
                resultEntity = getName(element)
                result = resultEntity == entity
                numResults += result
            if numResults == 0:
                forwardMisses += 1
                Input = input('A reverse lookup failed. Show context? y/n')
                if Input.lower() == 'y':
                    print('While mapping alias %s back to entity %s, the result was instead %s' % (alias, entity, resultEntity))
                Input = input('Exit test? y/n')
                if Input.lower() == 'y':
                    gate = True
                    break
            if gate:
                break
        bar.next()
        if gate:
            break
    bar.finish()
    print('%d look-ups failed.' % forwardMisses)
    print('The forward look-up dictionary had %d entries. %d were tested.' % (numEntities, entityCount))

    # reverseDic
    numAliases = len(reverseDic)

    # Progress bar
    if True:
        timeStamp = ts()
        print('\n[%s] Testing all reverse lookup dictionary contents' % str(timeStamp))
        count = numAliases
        bar = ChargingBar('', max = count)
        Tic = tic()

    # Test all reverse lookup dictionary contents
    reverseMisses = 0
    aliasCount = 0
    gate = False
    for alias, list in reverseDic.items():
        aliasCount += 1
        for element in list:
            entity = getName(element)
            reverseList = forwardDic[entity]
            numResults = 0
            for element in reverseList:
                resultAlias = getName(element)
                result = resultAlias == alias
                numResults += result
            if numResults == 0:
                reverseMisses += 1
                Input = input('A reverse lookup failed. Show context? y/n')
                if Input.lower() == 'y':
                    print('While mapping entity %s back to alias %s, the result was instead %s' % (entity, alias, resultAlias))
                Input = input('Exit test? y/n')
                if Input.lower() == 'y':
                    gate = True
                    break
            if gate:
                break
        bar.next()
        if gate:
            break
    bar.finish()
    print('%d reverse look-ups failed.' % reverseMisses)
    print('The reverse look-up dictionary has %d entries. %d were tested.' % (numAliases, aliasCount))

    if (forwardMisses + reverseMisses) == 0:
        print('The tested dictionaries are complete in both the forward and backward directions.')
    else:
        print('Warning. The dictionaries failed some forward or reverse lookups.')

# Test if values for s and m keys are same

def getCidSet(cpiDic):
    '''
    cpiDic --> cidSet

    Returns the CID number without prefix for all CID keys in the CPI Dictionary.

    Example:
    --------
    >>> cpiDic = {'CIDm53255435' : ['9606.ENSP00000256119', '9606.ENSP00000285379']}
    >>> getCidSet(cpiDic)
    {53255435}
    '''
    cidList = list(cpiDic.keys())
    cidSet = set()

    for cid in cidList:
        cid = cid[4:]
        cidSet.add(cid)

    return cidSet

# Get snapshot of cpiDic CIDs
def countCidTypes1(cpiDic):
    '''
    cpiDic --> counts of CIDm and CIDs values, but not unique and equivalent CID numbers. See getIsomerSetCounts.
    '''
    cpiKeys = cpiDic.keys()
    mCount = 0
    sCount = 0
    for cid in cpiKeys:
        if cid[:4] == 'CIDm':
            mCount += 1
        elif cid[:4] == 'CIDs':
            sCount += 1
    allCounts = mCount + sCount
    total = len(cpiKeys)
    missing = total - allCounts
    print('mCount: %d\nsCount: %d\ntotal:  %d\nmissing: %d' % (mCount, sCount, total, missing))
    return mCount, sCount, total, missing

def countCidTypes2(cpiDic):
    '''
    cpiDic --> counts unique

    Illustration
    ------------
    CID  |  CIDs    | CIDm
    -------------------------
    ..1  |  CIDs..1 | CIDm..1
    ..2  |  ---     | CIDm..2
    ..3  |  CIDs..3 | ---
    ..4  |  CIDs..4 | CIDm..4
    '''
    cpiKeys = cpiDic.keys()
    d = {}
    mSet = set()
    sSet = set()
    msSet = set()
    for cid in cidSet:
        mCheck = 'CIDm'+cid in cpiKeys
        sCheck = 'CIDs'+cid in cpiKeys
        if mCheck and sCheck:
            d[cid] = 'ms'
            msSet.add(cid)
        elif mCheck:
            d[cid] = 'm'
            mSet.add(cid)
        elif sCheck:
            d[cid] = 's'
            sSet.add(cid)
    lms = len(msSet)
    lm = len(mSet)
    ls = len(sSet)
    text = \
    '%d | Number of merged stereo-isomers (\'m\' CIDs)\n\
%d | Number of stereo isomers (\'s\' CIDs)\n\
%d | Number of \'m\' and \'s\' CIDs' % (lm, ls, lms)
    print(text)
    return lm, ls, lms, msSet, mSet, sSet

def countCidTypes(cpiDic):
    mCount, sCount, tCount, missing = countCidTypes1(cpiDic)
    lm, ls, lms, msSet, mSet, sSet = countCidTypes2(cpiDic)

    if False:
        print(mCount == lms+lm)
        print(sCount == lms+ls)
        print(tCount == 2*lms+lm+ls)

    a = mCount == lms+lm
    b = sCount == lms+ls
    c = tCount == 2*lms+lm+ls

    return a and b and c

def f():
    pass

    lm, ls, lms, msSet, mSet, sSet = countCidTypes2(cpiDic)

    for cid in msSet:


# Are m and s CIDs equivalent?
results = 0
for cid in cidSet:
    a = cpiDic['CIDm'+cid]
    b = cpiDic['CIDs'+cid]
    if a == b:
        results += 1
print(results)
