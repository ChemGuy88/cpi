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
