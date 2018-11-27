# Lab Journal for CPI

## To-do list
- [ ] Test `makeCpiDic`
- [ ] Test `makeChemSynsDic`
- [ ] Fix CID prefix issue

## Wishlist
- [ ] Should consider **differentiating DIR** to (DIR, outDIR, dataDIR), where
  a. DIR is the working directory containing cpi.py,
  b. outDIR is the directory where results are output (e.g., the pickled dictionaries), and
  c. dataDIR is the directory where primary data is located (e.g., the STRING and STITCH files)
- [ ] **Generalize `make` functions.**  
Some functions are exactly the same except for the core of the loop which does a different operation on each line of the links file. I should create a main function that takes another function to operate on each line. E.g.,:
  ```python
  for line in lines:
      output = function(line)
  ```
- [ ] **Remove usage of loadFile**  
Remove usage of `loadFile` from
  * makeProtSynsDic
  * makeLinksDic
- [ ] **Use `getQuickModeLimit` on applicable functions**  
  - [ ] makeCidList
  - [ ] makeProtSynsDic
  - [ ] makeLinksDic
  - [ ] loadFile  
- [ ] Iteratively import modules on `cpirun.py`. See this [link](https://www.tutorialspoint.com/Can-we-iteratively-import-python-modules-inside-a-for-loop).
- [ ] Display molecule images on a second column besides the tables from (#11/23/18)

## Journal Entries

* [Monday 11/26](#11/26/18)
* [Sunday 11/25](#11/25/18)
* [Saturday 11/24](#11/24/18)
* [Friday 11/23](#11/23/18)
* [Wednesday 11/21](#11/21/18)
* [Tuesday 11/20](#11/20/18)
* [Monday 11/19](#11/19/18)
* [Friday 11/9](#11/9/18)
* [Tuesday 11/6](#11/6/18)
* [Monday 10/29](#10/29/18)
* [Tuesday 10/23](#10/23/18)
* [Monday 10/22](#10/22/18)
* [Friday 10/19](#10/19/18)
* [Thursday 10/11](#10/11/18)
* [Tuesday 10/9](#10/9/18)
* [Friday 10/5](#10/5/18)
* [Monday 10/1](#10/1/18)

## [Completed To-do list](#CompletedToDoList)

### <a name="11/27/18"></a> Tuesday 11/27

Project put on hold in favor of a word-embedding project.

### <a name="11/26/18"></a> Monday 11/26

PubChem can have separate entries for left and right-handed enantiomers for some of their compounds, but not all. For example, [Thalidomide](https://www.acs.org/content/acs/en/molecule-of-the-week/archive/t/thalidomide.html) has three entries, one each for the left and right-handed enantiomers, and one without a specified stereochemistry. Another chemical, bromochlorofluoroethane, does not have such special treatment, and instead has only one entry with an undefined stereochemistry. I don't know if that's because there's no known synthetic path to isolate the enantiomers, i.e., the enantiomers don't exist. Consider also bromochlorofluoroiodomethane, which according to Wikipedia is a "hypothetical" chemical, and bromochlorofluoromethane. The former is listed in PubChem, but only with an undefined stereochemistry, while the latter has all three possible entries in PubChem, the two defined stereoisomers and the one undefined stereoisomer.

| Compound | Compound ID | Enantiomers Available on PubChem |
|----------|-------------|----------------------------------|
| Thalidomide     | [`00005426`](https://pubchem.ncbi.nlm.nih.gov/compound/5426#section=Top) | Yes |
| (S)-Thalidomide | [`00092142`](https://pubchem.ncbi.nlm.nih.gov/compound/92142#section=Top) | Yes |
| (R)-Thalidomide | [`00075792`](https://pubchem.ncbi.nlm.nih.gov/compound/75792#section=Top) | Yes |
| bromochlorofluoroethane | [`57450899`](https://pubchem.ncbi.nlm.nih.gov/compound/57450899#section=Top) | No |
| | | |
| Bromochlorofluoromethane | [`00079058`](https://pubchem.ncbi.nlm.nih.gov/compound/79058) | Yes |
| (S)-Bromochlorofluoromethane | [`57518772`](https://pubchem.ncbi.nlm.nih.gov/compound/57518772) | Yes |
| (R)-Bromochlorofluoromethane | [`57518771`](https://pubchem.ncbi.nlm.nih.gov/compound/57518771) | Yes |
| Bromo(Chloro)Fluoro(Iodo)Methane | [`57424940`](https://pubchem.ncbi.nlm.nih.gov/compound/57424940) | No |

Here's something else of interest. Searching the CID of all three Thalidomide entries in `cpidDic` gives the following results:

```Python
Thalidomide_x_m = cpiDic['CIDm00005426']
Thalidomide_x_s = cpiDic['CIDs00005426']
# Thalidomide_s_m = cpiDic['CIDm00092142'] # KeyError
Thalidomide_s_s = cpiDic['CIDs00092142']
# Thalidomide_r_m = cpiDic['CIDm00075792'] # KeyError
# Thalidomide_r_s = cpiDic['CIDs00075792'] # KeyError

Thalidomide_x_m == Thalidomide_x_s # True
Thalidomide_x_m == Thalidomide_s_s # False
```

After putting off the CID prefix issue I started working on the dictionaries again. While rewriting the `makeChemSynsDic` function it occurred to me that that function was overwriting results for the same CID. I wonder if it will be different after I finish rewriting it. I.e., the old cold was:

```python
chemSynsDic[cid] = listOfSynonyms
```

instead of

```python
try:
    chemSynsDic[cid].append(resultResult)
except KeyError:
    chemSynsDic[cid] = [resultResult]
```

Also updated `makeCpiDic` to return reverse dictionary.

Next, I need to test `makeCpiDic` and `makeChemSynsDic` to make sure the updates work.

### <a name="11/25/18"></a> Sunday 11/25

At first glance the bags have no apparent pattern.

```
bag: mseq
  mean:
  [1.4 1.2 0.2 0.1 0.1 0. ]
  min:
  [0 0 0 0 0 0]
  max:
  [4 4 1 1 1 0]
  Number of key errors: 1
  Bad keys: 05318656

bag: ms
  mean:
  [2.9 1.7 1.2 0.2 0.2 0. ]
  min:
  [0 0 0 0 0 0]
  max:
  [11  6  5  1  1  0]
  Number of key errors: 7
  Bad keys: 09861986, 00065153, 00456846, 01819399, 03109913, 06918304, 00017793

bag: m
  mean:
  [4.9 4.5 0.4 0.3 0.3 0. ]
  min:
  [0 0 0 0 0 0]
  max:
  [17 17  2  1  1  0]
  Number of key errors: 3
  Bad keys: 00090173, 00102175, 05284237

bag: s
  mean:
  [5.6 3.8 1.8 0.  0.  0. ]
  min:
  [2 0 0 0 0 0]
  max:
  [20 20  8  0  0  0]
  Number of key errors: 5
  Bad keys: 00044568, 00130001, 05351908, 00007528, 00000207
```

### <a name="11/24/18"></a> Saturday 11/24

Edited `downloadCidSyns` to also scrape Computed Physical Properties from PUG REST. Edited `makeChemSynsDic` to parse the corresponding scraped information. I just have to analyze it to see if a pattern emerges that corresponds to the m and s labels.

### <a name="11/23/18"></a> Friday 11/23

PubChem has useful information about a compound's stereochemistry under "4 Chemical and Physical Properties: 4.1 Computed Properties". I can thus programmatically use this information to elucidate some pattern in which the "m" and "s" labels are assigned.

#### Tables of Groupings

What follows is a series of tables. Each table has 5 compounds from each possible CID label combination described in the table from [Wednesday 11/19](#11/19/18). "None" under "Stereochemistry" means that it has no chiral centers and so the molecule has no stereoisomers.

###### Numbers that are both m and s-types **AND** point to the same set of proteins (`mseq`)

| CID Number | Stereochemistry | Note |
|------------|-----------------|------|
| 00128493 | None | |
| 04348756 | None | |
| 25227472 | None | |
| 46230909 | None | |
| 66759619 | None | ||

<figure>
    <div style="text-align:center"><img width="200" height="200" align="middle" src='/labJournalDiagrams/mseq/00128493.png' />
    <font size="2" align="center"><figcaption> Figure 2: 00128493 </figcaption></font></div>
</figure>
<figure>
    <div style="text-align:center"><img width="200" height="200" align="middle" src='/labJournalDiagrams/mseq/04348756.png' />
    <font size="2" align="center"><figcaption> Figure 2: 04348756 </figcaption></font></div>
</figure>
<figure>
    <div style="text-align:center"><img width="200" height="200" align="middle" src='/labJournalDiagrams/mseq/25227472.png' />
    <font size="2" align="center"><figcaption> Figure 2: 25227472 </figcaption></font></div>
</figure>
<figure>
    <div style="text-align:center"><img width="200" height="200" align="middle" src='/labJournalDiagrams/mseq/46230909.png' />
    <font size="2" align="center"><figcaption> Figure 2: 46230909 </figcaption></font></div>
</figure>
<figure>
    <div style="text-align:center"><img width="200" height="200" align="middle" src='/labJournalDiagrams/mseq/66759619.png' />
    <font size="2" align="center"><figcaption> Figure 2: 66759619 </figcaption></font></div>
</figure>  

###### Numbers that are both m and s-types (`msbag`) **BUT** do not point to the same set of proteins

| CID Number | Stereochemistry | Note |
|----------|-------------------|------|
| 00073410 | Yes | "Conformer generation is disallowed since too many atoms" |
| 09801471 | Enantiomeric | |
| 11611496 | Yes | |
| 16129736 | Yes | "Conformer generation is disallowed since too many atoms" |
| 44273465 | Yes | ||


###### m-type CID numbers (`mbag`)

| CID Number | Stereochemistry | Note |
|----------|-------------------|------|
| 00002813 | None | |
| 00191957 | None | |
| 00557441 | Unknown | "Conformer generation is disallowed since too many undefined stereo centers" |
| 10162715 | Yes | Trans and a chiral center |
| 24951130 | Yes | ||

###### s-type CID numbers (`sbag`)

| CID Number | Stereochemistry | Note |
|----------|-------------------|------|
| 09576097 | Yes | |
| 44327287 | Enantiomeric | |
| 44370954 | Yes | |
| 54036959 | Yes | |
| 90832670 | Yes | ||

### <a name="11/21/18"></a> Wednesday 11/21

Fixed `getCidBags`.

### <a name="11/20/18"></a> Tuesday 11/20

Apparently my `getCidBags` function didn't work. The intersection of `mbag` and `sbag` has 15,000 numbers, or about 13% of their union (112,988).

### <a name="11/19/18"></a> Monday 11/19

I think the CID numbers that have equivalent m and s results (i.e., point to the same set of proteins), are compounds without stereoisomers, e.g., ethanol, gold, or NaCl. But then should they even have an 'm' number? I'm closer to this answer. Using the function `getCidBags` I am able to compare the protein sets pointed to by each m and s-type CID number. So far I only have count results. Next I will dive into the CID sets that have now been sorted. Below are the count results.

```
Table illustrating the possible values and combinations of m and s-type CID numbers:
 _______________________________________________________________
       m:        x |     None |        x |     None |        x |
       s:        x |     None |     None |        x |        y |

 Table with results of counting each of the above cases
 _______________________________________________________________
     var:     mseq |      ms0 |       m1 |       s1 |      ms1 |
   count:   329462 |        0 |    29697 |    68164 |    15127 |
```

The next step is to look at the Python sets returned by `getCidBags`, `mbag`, `sbag`, and `msbag` and see if there's any pattern to how the CID numbers have been sorted. I will have to look up the numbers manually in PUG REST.


### <a name="11/9/18"></a> Friday 11/9

The next thing to do is

1. compare `cpiDic[CIDs] == cpidDic[CIDm]` for all CID.
2. confirm that interaction for “m”-type chemical IDs implies interaction for all “s”-type chemical IDs.

>However, the issue is not completely resolved because (1) there are instances in which some “s” and “m” pairs do not point to the exact same set. (2) Some chemicals IDs are only “m” and some are only “s”. Also, (3) I have to confirm that interaction for “m”-type chemical IDs implies interaction for all “s”-type chemical IDs.

### <a name="11/6/18"></a> Tuesday 11/6

I had to clear my mind and start from scratch. The following are the core functions I need.
```Python
def makeCpiDic():
    '''
    9606.protein_chemical.links.v5.0.tsv --> cpiDic

    e.g.:
    key             value
    CIDm91758680    9606.ENSP00000257254
    '''
    pass

def getCidList():
    '''
    cpiDic --> cidList

    e.g.:
    ['CIDm53255435',
     'CIDs45273760',
     'CIDs11472526',
     'CIDs44423364',
     'CIDm71718190']
    '''
    pass

def makeChemSynsDic():
    '''
    same as previous version, but with prefixes
    '''
    pass

def makeProtSynsDic():
    '''
    same as previous version, but without loadFile
    '''
    pass
```

Fixing the prefixes issue is the main problem, so I kept working on that. I was just playing around with the dictionary to see what to do. Below is a summary of results.

```Python
mCount, sCount, tCount, missing = countCidTypes1(cpiDic)
# mCount: 412753
# sCount: 374286
# total:  787039
# missing: 0

lm, ls, lms = countCidTypes2(cpiDic)
# 29697 | Number of CIDs that appear just as merged stereo-isomers ('m' CIDs)
# 68164 | Number of CIDs that appear just as stereo isomers ('s' CIDs)
# 344589 | Number of CIDs that appear as both 's' and 'm'.

mCount == lms+lm
# True
sCount == lms+ls
# True
tCount == lms+lm+ls
# True
```

My understanding of the problem, I think, will be aided by elaborating on the meaning of the following visualization:
```
CID  |  CIDs    | CIDm
-------------------------
..1  |  CIDs..1 | CIDm..1
..2  |  ---     | CIDm..2
..3  |  CIDs..3 | ---
..4  |  CIDs..4 | CIDm..4
```

Besides that I also fleshed out the docstring to `cpi.py` and `cpirun.py`. I also added a new module, `cpiaux.py`. I moved all debugging functions to the latter.

### <a name="10/29/18"></a> Monday 10/29

I started by trying to see how important the CID prefixes are. In the set of all CIDs, about half are stereo-specific, and the other half are not. Until I learn more it seems they cannot be merged, because although they point to some common proteins, there are differences.

While working on the CID prefix issue I realized I can make the CPI dic by altering the `makeCidList` function. So I can put that new code under another function name, because I still need to make cidList to get all the aliases under PUG REST. I've put this new `make` function in `cpirun.py` under `makeCpiDic`. I have also saved the `cidList` and `cpiDic` objects as pickles in `\9606Pickles`. Because they took 12 and 21 minutes to generate, respectively.

### <a name="10/23/18"></a> Tuesday 10/23

Fixed the dictionary checkers by adding a new function `getName`. It is used for extracting `alias` and `entity` values.

Removed `quickMode` and `quickModeLimit` from `makeChemSynsDic` and `downloadCidSyns`, because `makeChemSynsDic` is a function of `cidLists` which is just a `list`, and it's easier to manipulate that list than to implement a `quickMode` for the function. However, if `makeChemSynsDic` at some point takes in a file that calls `makeCidList`, it will need values for `quickMode` and `quickModeLimit`.

Still need to:

1. fix the CID prefixes
2. check the dictionary makers again.
3. Test `getQuickModeLimit` on all functions that use it.

### <a name="10/22/18"></a> Monday 10/22

Continued working on fixing the `makeChemSynsDic` function. There are two problems.

1. The dictionary keys are in the format `123456789` instead of `CIDm123456789`/`CIDs123456789` format. This is because the synonyms are fetched from NIH via the PUG REST API. The `CIDm`/`CIDs` prefixes are not allowed in the API. So I have to restore the prefixes after the synonyms are fetched.
2. The dictionary checker doesn't work for `chemSynsDic` because the elements are string singletons, not tuples (as is the case in `protSynsDic`, for which I had originally designed the dictionary-checker function).

One possible solution for (2) is:

```Python
for element in aliases:
  alias = getAlias(element)
  pass
```

Another possible solution for the second problem is to replace

```Python
for alias, source in aliases:
  pass
```

with

```Python
for element in aliases:
  alias, source = [*element][0]
  pass
```

or

```Python
for element in aliases:
  if isinstance(element, string):
    alias = element
  elif isinstance(element, (tuple, list)):
    alias = element[0]
  pass
```

### <a name="10/19/18"></a> Friday 10/19

Here's an idea on how to build my corpus, based on [previous notes](#10/5/18).

1. Get list of all confirmed cpi
2. For each c-p pair, put sentences that mention the pair in a set, `cpi`
3. Get list of random c-p pairs that are not confirmed to interact.
4. For each c-p pair, put sentences that mention the pair in a set, `noncpi`
5. Operate on both sets
  1. The problem is reduced to a clustering machine learning problem.
  2. ...

To illustrate:  

![venDiagram](/labJournalDiagrams/venDiagram.png)


### <a name="10/11/18"></a> Thursday 10/11

Improving `makeChemSynsDic` function. Modified flowchart to match. I'm debating whether `makeChemSynsDic` can intake three files, and make `cidLists` from each one. I've made enough progress to splice the correction into the `cpi` module from my workspace.

Also, I think I realized why my method for counting lines in a file was not giving consistent results. I used `(f.readline for line in f)`. I suspect the reason this halves the number of lines is because `f.readline` happens after `for line in f`, so you're essentially reading every other line. Indeed, this is confirmed when I run the following code on the file `9606.protein_chemical.links.detailed.v5.0.tsv`.

```Python
In [112]: with open(fpath) as f:
     ...:     gen = (f.readline() for line in f)
     ...:     for _ in range(5):
     ...:         print(next(gen))
```

which yields

```
CIDm91758680	9606.ENSP00000257254	0	0	0	278	279

CIDm91758408	9606.ENSP00000006777	0	0	0	225	225

CIDm91758408	9606.ENSP00000216085	0	0	0	225	225

CIDm91758408	9606.ENSP00000235329	0	0	0	162	162

CIDm91758408	9606.ENSP00000262366	0	0	0	169	169
```

However, the first few lines (using bash command `more`), are actually

```
chemical        protein experimental    prediction      database        textmining      combined_score
CIDm91758680    9606.ENSP00000257254    0       0       0       278     279
CIDm91758680    9606.ENSP00000302120    0       0       0       154     154
CIDm91758408    9606.ENSP00000006777    0       0       0       225     225
CIDm91758408    9606.ENSP00000056217    0       0       0       178     178
CIDm91758408    9606.ENSP00000216085    0       0       0       225     225
CIDm91758408    9606.ENSP00000221740    0       0       0       151     151
CIDm91758408    9606.ENSP00000235329    0       0       0       162     162
CIDm91758408    9606.ENSP00000259791    0       0       0       194     194
CIDm91758408    9606.ENSP00000262366    0       0       0       169     169
```

If I count the number of newline characters `('\n')` in the same file, I get 21,773,492 (a). If I count the number of iterations it takes to do f.readline(), I get 10,886,746 (b). The ratio a/b is exactly 2, which again confirms my suspicion that the code was reading every other line.

### <a name="10/9/18"></a> Tuesday 10/9

Tested protein dictionaries (abridged versions) forward and backwards. They passed. This means that the `makeProtSynsDic` function should create a working full dictionary.

### <a name="10/5/18"></a> Friday 10/5

Notes from my presentation to the lab group.

1. Use Xin's protein-protein interaction (ppi) model for `protein_A`-`protein_B` interactions, but instead of using `protein_A`, use `chemical_C`. This should work because the syntax between ppi and cpi corpus should be similar. What the neural network is doing looking for syntactical patterns in the training set. The way Xin has written his code the training operates over a set of text with the `protein` value substituted with a place-holder. The actual values of `protein_A` don't matter; it's the words revolving around the `protein` variables.
2. Terms related to this kind of research:
  1. Transfer learning
  2. Distance learning
  3. Multiple instance learning
3. There is a python module that parses sentences into grammatical components. It might be called NSTK. Related terms: "Parse dependency approach," "dependency tree".
4. Based on Dr. Zhang's explanation of his previous work in NLP (the Harvard paper), I'm under the impression that the neural network is finding the intersection of a series of sets. This intersection represents the sentences that have a syntax that implies chemical-protein interaction.


### <a name="10/1/18"></a> Monday 10/1

Finished testing the protein synonyms dictionary with the CPI dictionary. All dictionaries were limited to just species 9606 (human).

```
======================================================================
CPI dictionary has 0 unique proteins.
protein synonyms dictionary has 3105 unique proteins.
They share 17352 (84.82%) proteins.
======================================================================
```

It is obvious from these results that all the interaction proteins (from species 9606) are in the synonyms dictionary (from species 9606). But not all proteins are said to have interactions. This raises the question of whether the dictionary is complete or was made incorrectly.

## <a name="CompletedToDoList"></a> Completed To-do list
- [x] **Place `TicSum` in `downloadCidSyns` under for-loop**
- [ ] Fix CID prefix issue
  - [x] See if the CIDs from `countCidTypes2` that are in the `ms` bag point to the same protein sets.
