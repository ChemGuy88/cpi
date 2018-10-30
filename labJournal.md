# Lab Journal for CPI

## To-do list

1. Should consider **differentiating DIR** to (DIR, outDIR, dataDIR), where
  a. DIR is the working directory containing cpi.py,
  b. outDIR is the directory where results are output (e.g., the pickled dictionaries), and
  c. dataDIR is the directory where primary data is located (e.g., the STRING and STITCH files)
2. **Generalize `make` functions.**  
Some functions are exactly the same except for the core of the loop which does a different operation on each line of the links file. I should create a main function that takes another function to operate on each line. E.g.,:

  ```python
  for line in lines:
      output = function(line)
  ```
3. **Remove usage of loadFile**  
Remove usage of `loadFile` from
  * makeProtSynsDic
  * makeLinksDic

4. **Use `getQuickModeLimit` on applicable functions**  
  * loadFile  
  * makeCidList
  * makeProtSynsDic
  * makeLinksDic

## Journal Entries

* [Monday 10/29](#10/29/18)
* [Tuesday 10/23](#10/23/18)
* [Monday 10/22](#10/22/18)
* [Friday 10/19](#10/19/18)
* [Thursday 10/11](#10/11/18)
* [Tuesday 10/9](#10/9/18)
* [Friday 10/5](#10/5/18)
* [Monday 10/1](#10/1/18)

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

It is obvious from these results that all the interaction proteins (from species 9606) are in the synonyms dictionary (from species 9606). But not all proteins are said to have interactions. This raises the
