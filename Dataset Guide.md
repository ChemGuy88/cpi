Guide to the Chemical-Protein Interaction Datasets
==================================================

## Note

**All files and sites are STRING version 10, unless otherwise noted.**

## Purpose

This document is just me putting my thoughts down to paper to better understand the datasets. I hope this helps someone else, too.

## Contents

1. [The Chemical-Protein (Inter)actions File](#CPLinks)
   1. [Chemical ID Prefix](#CIDPrefix)
   2. [Detail 2](#Detail2)
1. [The Protein Dictionary](#ProteinDic)
2. [The Species Dictionary](#SpeciesDic)
3. [The Chemicals Dictionary](#ChemicalsDic)
4. [Other Files](#OtherFiles)
5. [Comments](#Comments)

### <a name="CPLinks"></a> 1. The Chemical-Protein (Inter)actions File, or The Only File You Need for Now

The Chemical-Protein pairs and their interactions are listed in `9606.actions.v5.0.tsv`. Note that this particular file name is prepended with `9606`, which denotes links for just one species, Humans. The complete dataset is massive, at around 60+ Gigabytes. Also note that this file contains bidirectional graph data. So each link is listed twice, for the forward and backward direction.

Preview [(1)](#noteFileLegibility):

```
item_id_a             item_id_b             mode        action  a_is_acting score
9606.ENSP00000170630  CIDm00010461          expression          f           150
CIDm00010461          9606.ENSP00000170630  expression          t           150
9606.ENSP00000353915  CIDs23627457          binding             f           191
CIDs23627457          9606.ENSP00000353915  binding             f           191
```

#### 1. <a name="CIDPrefix"></a> Chemical ID prefix

   The chemical IDs are the values that begin with `CID`[(2)](#noteCIDvsSID). According to the [STITCH README](http://stitch.embl.de/download/README):

>CIDs / CID0... - this is a stereo-specific compound, and the suffix is the
PubChem compound id.
>
>CIDm / CID1... - this is a "flat" compound, i.e. with merged stereo-isomers
The suffix (without the leading "1") is the PubChem compound id.
>
>Note that the download files contain the prefixes CIDs/CIDm, while the API
still returns CID0/CID1.)

#### 2. <a name="Detail2"></a> Detail 2

Fill me up, butter cup.

### 2. The Chemical-Protein **Interactions** Dictionary, or lack thereof

The Chemical-Protein interaction list from my database (STITCH) has only 7 types of interactions (listed below). That means I won't need a dictionary for them. For example, when the neural net reads an abstract, the interaction term "binding" is just regular English, so it can use word embedding for it. None of the interaction terms seem to be technical terms. It's not like "acetaminophen", which needs to be translated (with a dictionary), to Tylenol or something else. The interaction terms extracted are a set of 7
  1. `activation`
  2. `binding`
  3. `catalysis`
  4. `expression`
  5. `inhibition`
  6. `reaction`
  7. `pred_bind`

### <a name="ProteinDic"></a> 1. The Protein Dictionary

The `protein.aliases.v[...].txt` files seem to be protein dictionaries. They have the generic name on the leftmost column, the alias in the seonc column, and the source of that alias in the third column. `[...]` is the file version number. Currently I'm using version `10.5`.

### <a name="SpeciesDic"></a>2. The Species Dictionary

A file or other object that links STITCH-format **species** identifiers to any of its synonyms used in the literature.

### <a name="ChemicalsDic"></a>3. The Chemicals Dictionary

A file or other object that links STITCH-format **chemicals** identifiers to any of its synonyms used in the literature. I have made my own using the NIH's PUG REST service.

### 4. <a name="OtherFiles"></a>Other Files

#### 1. Other File 1

Text.

#### 2. Other File 2

Text

### 5. <a name="Comments"></a>Comments

1. <a name="noteCIDvsSID"></a>The prefix `CID` actually stands for "Compound ID", and should be distinguinshed from another chemical ID prefix (not in this dataset), `SID`, which stands for "Substance ID". These are PubChem conventions.
1. <a name="noteFileLegibility"></a> The file is not spaced as shown. The file is tab-delimited or space delimited. The spacing was added for better legibility.
