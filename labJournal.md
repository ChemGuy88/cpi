# Lab Journal for CPI

## To-do list

1. Make a general-purpose function that checks forward and reverse dictionaries. Use the workspace code from `cpirun.py`.
1. Should consider differentiating DIR to (DIR, outDIR, dataDIR), where
  a. DIR is the working directory containing cpi.py,
  b. outDIR is the directory where results are output (e.g., the pickled dictionaries), and
  c. dataDIR is the directory where primary data is located (e.g., the STRING and STITCH files)
2. **Generalize `make` functions.**  
Some functions are exactly the same except for the core of the loop which does a different operation on each line of the links file. I should create a main function that takes another function to operate on each line. E.g.,:

  ```python
  for line in lines:
      output = function(line)
  ```
4. **Check lower bound for `quickModeLimit` in...**  
  ...
3. **Add fail-safe for `quickModeLimit`**.  
  `quickModeLimit` is set to a specific value for some functions. A fail-safe should be added in case the `quickModeLimit` is greater than the number of iterable items.

## Journal Entries

* [Thursday 10/11](#10/11/18)
* [Tuesday 10/9](#10/9/18)
* [Friday 10/5](#10/5/18)
* [Monday 10/1](#10/1/18)

### <a name="10/19/18"></a> Monday 10/22

Continued working on fixing the `makeChemSynsDic` function. There are two problems.

1. The dictionary keys are in the format `123456789` instead of `CIDm123456789`/`CIDs123456789` format. This is because the synonyms are fetched from NIH via the PUG REST API. The `CIDm`/`CIDs` prefixes are not allowed in the API. So I have to restore the prefixes after the synonyms are fetched.
2. The dictionary checked doesn't work for `chemSynsDic` because the elements are string singletons, not tuples (as is the case in `protSynsDic`, for which I had originally designed the dictionary-checker function).

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
  if len(element) == 0:
      alias = element
  else:
      alias = [*element][0]
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

![venDiagram](/Users/Herman/Documents/jzhang/cpiDiagrams/venDiagram.png)


### <a name="10/11/18"></a> Thursday 10/11

Improving `makeChemSynsDic` function. Modified flowchart to match. I'm debating whether `makeChemSynsDic` can intake three files, and make `cidLists` from each one. I've made enough progress to splice the correction into the `cpi` module from my workspace.
Also, I think I realized why my method for counting lines in a file was not giving consistent results. I used `(f.readline for line in f)`. I suspect the reason this halves the number of lines is because `f.readline` happens after `for line in f`, so you're essentially reading every other line. Indeed, this is confirmed when I run
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

### <a name="10/9/18"></a> Tuesday 10/9

Tested protein dictionaries (abridged versions) forward and backwards. They passed. This means that the `makeDic` function should create a working full dictionary.

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
