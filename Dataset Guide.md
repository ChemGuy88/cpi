Guide to the Chemical-Protein Interaction Datasets
==================================================

## Note

**The STRING site cited is version 7, and version 10 is available.**

## Purpose

This document is just me putting my thoughts down to paper to better understand the datasets. I hope this helps someone else, too.

## Contents

1. [The Protein Dictionary](#ProteinDic)
2. [The Species Dictionary](#SpeciesDic)
3. [The Chemicals Dictionary](#ChemicalsDic)
4. [Other Files](#OtherFiles)
5. [Comments](#Comments)

### <a name="ProteinDic"></a> 1. The protein Dictionary

According to the Readme from STITCH (http://stitch1.embl.de/download/README), the protein dictionary can be found in another database, STRING. The file is `protein.aliases.v7.1.txt`. Its first two lines are[(1)](#noteFileLegibility):

```
##  species_ncbi_taxon_id ##  protein_id ## alias ##  source ##
    117                       RB2638        1,4-alpha D-glucan:1,4-alpha D-glucan 6-glucosyl- transferase	paralign_SWISSPROT_DE
```

However, when queried on RefSeq, this tag does not return the same aliases, but rather simply "glycogen branching enzyme".

Looking around STRING we see that they have a map of their [datasets](http://string71.embl.de/newstring_download/database.schema.v7.1.pdf). In it they explain that the protein tag comes from third party databases: GenomeReviews, RefSeq, or Ensembl. However, GenomeReviews is now defunct and forwards to Ensembl.

As an example, consider the CPI pair (in Python dictionary format):

```python
'CID000007215' : ('117.RB7803', '229')
```

`117` is the species, and `RB7803` is the protein.

The protein can be found in [RefSeq](https://www.ncbi.nlm.nih.gov). It actually points to a DNA sequence, but it has a protein associated with it (uroporphyrinogen III synthase, uroporhyrinogen decarboxylase)

For another CPI pair,
```python
'CID006918553' : ('7955.ENSDARP00000041069', '278')
```

`ENSDARP00000041069` points to a DNA sequence on [Ensembl](http://useast.ensembl.org). This only points to a DNA sequence, but the search result yields "ENSDARG00000032319", which is associated with the protein "mest". `ENSDARP00000041069`, however, does return results for a DNA sequence that is associated with "mest". I do not understand why this protein tag does not link directly to its associated gene.

### <a name="SpeciesDic"></a>2. The Species Dictionary

The species dictionary is located in file the `species.v7.1.txt`.

The first two lines are:

```
##  taxon_id  STRING_type STRING_name_compact     official_name_NCBI
    117       core        Rhodopirellula baltica  Rhodopirellula baltica (strain 1)
```

### <a name="ChemicalsDic"></a>3. The Chemicals Dictionary

The chemicals dictionary is in `chemicals.v1.10.tsv`. The first two lines are:

```
chemical      name                molecular_weight    SMILES_string
CID000000001  acetyl-L-carnitine  203.236             CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C
```

### 4. <a name="OtherFiles"></a>Other Files

   1. The protein-chemical links list, `protein_chemical.links.v1.0.tsv`:

```
chemical      protein                 combined_score
CID011989247  9031.ENSGALP00000004303 268
```

  The newer version for human links (`9606.protein_chemical.links.v5.0.tsv`), has the following preview:

```
chemical	protein	combined_score
CIDm91758680	9606.ENSP00000257254	279
CIDm91758680	9606.ENSP00000302120	154
```

  According to the [STITCH README](http://stitch.embl.de/download/README):

>CIDs / CID0... - this is a stereo-specific compound, and the suffix is the
PubChem compound id.
>
>CIDm / CID1... - this is a "flat" compound, i.e. with merged stereo-isomers
The suffix (without the leading "1") is the PubChem compound id.
>
>Note that the download files contain the prefixes CIDs/CIDm, while the API
still returns CID0/CID1.)

2. The link type list, `9606.actions.v5.0.tsv`. Note that this particular file name is prepended with `9606`, which denotes links for just one species, Humans. The complete dataset is massive, at around 60+ Gigabytes. Also note that this file contains bidirectional graph data. So each link is listed twice, for the forward and backward direction.

```
item_id_a             item_id_b             mode        action  a_is_acting score
9606.ENSP00000170630  CIDm00010461          expression          f           150
CIDm00010461          9606.ENSP00000170630  expression          t           150
9606.ENSP00000353915  CIDs23627457          binding             f           191
CIDs23627457          9606.ENSP00000353915  binding             f           191
```

### 5. <a name="Comments"></a>Comments

1. <a name="noteFileLegibility"></a> The file is not spaced as shown. The file is tab-delimited or space delimited. The spacing was added for better legibility.
