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

### Thursday 10/11

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

### Tuesday 10/9

Tested protein dictionaries (abridged versions) forward and backwards. They passed. This means that the `makeDic` function should create a working full dictionary.

### Monday 10/1

Finished testing the protein synonyms dictionary with the CPI dictionary. All dictionaries were limited to just species 9606 (human).

```
======================================================================
CPI dictionary has 0 unique proteins.
protein synonyms dictionary has 3105 unique proteins.
They share 17352 (84.82%) proteins.
======================================================================
```

It is obvious from these results that all the interaction proteins (from species 9606) are in the synonyms dictionary (from species 9606). But not all proteins are said to have interactions. This raises the
