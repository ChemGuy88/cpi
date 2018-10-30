Guide to the Chemical-Protein Interaction Datasets
==================================================

## Note

**All files and sites are from STITCH version 5, unless otherwise noted.**

## Purpose

This document is just me putting my thoughts down to paper to better understand the datasets. I hope this helps someone else, too.

## Contents

1. [Summary of Files and Contents](#fileTypeSummary)
1. [The Chemical-Protein Interactions File](#cpiDic)
    * [Chemical ID Prefix](#CIDPrefix)
2. [The Chemical-Protein Interaction _Types_ File](#cpiTypesDic)
3. [The Protein Dictionary](#ProteinDic)
4. [The Species Dictionary](#SpeciesDic)
5. [The Chemicals Dictionary](#ChemicalsDic)
6. [Other Files](#OtherFiles)
7. [Comments](#Comments)

### <a name="fileTypeSummary"></a> 1. Summary of Files and Contents

`9606.actions.v5.0.tsv`
```
item_id_a             item_id_b             mode        action  a_is_acting score
9606.ENSP00000170630  CIDm00010461          expression          f           150
CIDm00010461          9606.ENSP00000170630  expression          t           150
9606.ENSP00000353915  CIDs23627457          binding             f           191
CIDs23627457          9606.ENSP00000353915  binding             f           191
```
`9606.protein_chemical.links.v5.0.tsv`
```
chemical        protein                 combined_score
CIDm91758680    9606.ENSP00000257254    279
CIDm91758680    9606.ENSP00000302120    154
CIDm91758408    9606.ENSP00000006777    225
CIDm91758408    9606.ENSP00000056217    178
CIDm91758408    9606.ENSP00000216085    225
CIDm91758408    9606.ENSP00000221740    151
```
`9606.protein_chemical.links.detailed.v5.0.tsv`
```
chemical        protein experimental    prediction      database        textmining      combined_score
CIDm91758680    9606.ENSP00000257254    0               0               0               278     279
CIDm91758680    9606.ENSP00000302120    0               0               0               154     154
CIDm91758408    9606.ENSP00000006777    0               0               0               225     225
CIDm91758408    9606.ENSP00000056217    0               0               0               178     178
CIDm91758408    9606.ENSP00000216085    0               0               0               225     225
CIDm91758408    9606.ENSP00000221740    0               0               0               151     151
```
`9606.protein.aliases.v10.5.txt`
```
## string_protein_id    ## alias                                                ## source ##
9606.ENSP00000218548    #945                                                    BLAST_KEGG_NAME
9606.ENSP00000282841    (2E,6E)-farnesyl diphosphate synthase                   BLAST_UniProt_DE Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DE Ensembl_UniProt_DE
9606.ENSP00000349078    (2E,6E)-farnesyl diphosphate synthase                   BLAST_UniProt_DE Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DE Ensembl_UniProt_DE
9606.ENSP00000215574    (E3-independent) E2 ubiquitin-conjugating enzyme R1     BLAST_UniProt_DE Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DE Ensembl_UniProt_DE
9606.ENSP00000303709    (E3-independent) E2 ubiquitin-conjugating enzyme E1     BLAST_UniProt_DE Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DE Ensembl_UniProt_DE
```
`chemicals.v1.0.tsv`
```
chemical        name                                molecular_weight   SMILES_string
CID000000001    acetyl-L-carnitine                  203.236            CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C
CID000000003    2,3-dihydro-2,3-dihydroxybenzoate   156.136            C1=CC(C(C(=C1)C(=O)O)O)O
CID000000004    1-amino-2-propanol                  75.1097            CC(CN)O
CID000000005    3-amino-2-oxopropyl phosphate       169.073            C(C(=O)COP(=O)(O)O)N
CID000000006    DNCB                                202.552            C1=CC(=C(C=C1[N+](=O)[O-])[N+](=O)[O-])Cl
CID000000007    9-ethyladenine                      163.18             CCN1C=NC2=C1N=CN=C2N
```


### <a name="cpiDic"></a> 2. The Chemical-Protein Interactions File

The Chemical-Protein pairs and their interactions are listed in `9606.protein_chemical.links.v5.0.tsv`. Note that this particular file name is prepended with `9606`, which denotes links for just one species, Humans. The complete dataset is massive, at around 30.9 Gigabytes.

The following is a preview of the STRING file [(1)](#noteFileLegibility):

```
chemical        protein                 combined_score
CIDm91758680    9606.ENSP00000257254    279
CIDm91758680    9606.ENSP00000302120    154
CIDm91758408    9606.ENSP00000006777    225
CIDm91758408    9606.ENSP00000056217    178
CIDm91758408    9606.ENSP00000216085    225
CIDm91758408    9606.ENSP00000221740    151
```

##### <a name="CIDPrefix"></a> Chemical ID prefix

The chemical IDs are the values that begin with `CID`[(2)](#noteCIDvsSID). According to the [STITCH README](http://stitch.embl.de/download/README):

>CIDs / CID0... - this is a stereo-specific compound, and the suffix is the
PubChem compound id.
>
>CIDm / CID1... - this is a "flat" compound, i.e. with merged stereo-isomers
The suffix (without the leading "1") is the PubChem compound id.
>
>Note that the download files contain the prefixes CIDs/CIDm, while the API
still returns CID0/CID1.)

### <a name="cpiTypesDic"></a> 3. The Chemical-Protein Interaction _Types_ File.

According to the [STITCH README](http://stitch.embl.de/download/README), the interaction types are a subset of the total number of interactions. In other words, some CP pairs are known to interact, but their interaction type is not specified. STITCH stores these values in `9606.actions.v5.0.tsv`.

The interaction terms extracted are a set of 7
  1. `activation`
  2. `binding`
  3. `catalysis`
  4. `expression`
  5. `inhibition`
  6. `reaction`
  7. `pred_bind`

Note that this file contains bidirectional graph data. So each link is listed twice, for the forward and backward direction.

The following is a preview of the STRING file [(1)](#noteFileLegibility):

```
item_id_a             item_id_b             mode        action  a_is_acting score
9606.ENSP00000170630  CIDm00010461          expression          f           150
CIDm00010461          9606.ENSP00000170630  expression          t           150
9606.ENSP00000353915  CIDs23627457          binding             f           191
CIDs23627457          9606.ENSP00000353915  binding             f           191
```

### <a name="ProteinDic"></a> 4. The Protein Dictionary

A file or other object that links STRING-format **protein** identifiers to any of its synonyms used in the literature. STRING stores its protein dictionary in `protein.aliases.v10.5.txt`. They have the generic name on the leftmost column, the alias in the second column, and the source of that alias in the third column. Note this is the only file not from STITCH.

The following is a preview of the STRING file [(1)](#noteFileLegibility):

```
## string_protein_id    ## alias                                                ## source ##
9606.ENSP00000218548    #945                                                    BLAST_KEGG_NAME
9606.ENSP00000282841    (2E,6E)-farnesyl diphosphate synthase                   BLAST_UniProt_DE Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DE Ensembl_UniProt_DE
9606.ENSP00000349078    (2E,6E)-farnesyl diphosphate synthase                   BLAST_UniProt_DE Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DE Ensembl_UniProt_DE
9606.ENSP00000215574    (E3-independent) E2 ubiquitin-conjugating enzyme R1     BLAST_UniProt_DE Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DE Ensembl_UniProt_DE
9606.ENSP00000303709    (E3-independent) E2 ubiquitin-conjugating enzyme E1     BLAST_UniProt_DE Ensembl_HGNC_UniProt_ID(supplied_by_UniProt)_DE Ensembl_UniProt_DE
```

### <a name="SpeciesDic"></a> 5. The Species Dictionary

A file or other object that links STITCH-format **species** identifiers to any of its synonyms used in the literature. Currently there is no need to implement this.

### <a name="ChemicalsDic"></a> 6. The Chemicals Dictionary

A file or other object that links STITCH-format **chemicals** identifiers to any of its synonyms used in the literature. STITCH stores its chemical dictionary in `chemicals.v5.0.tsv.` I have made my own using the NIH's PUG REST service. STITCH offers a dictionary, although with only one alias. The PUG REST version has more aliases.

The following is a preview of the STITCH file [(1)](#noteFileLegibility):

```
chemical        name                                molecular_weight   SMILES_string
CID000000001    acetyl-L-carnitine                  203.236            CC(=O)OC(CC(=O)[O-])C[N+](C)(C)C
CID000000003    2,3-dihydro-2,3-dihydroxybenzoate   156.136            C1=CC(C(C(=C1)C(=O)O)O)O
CID000000004    1-amino-2-propanol                  75.1097            CC(CN)O
CID000000005    3-amino-2-oxopropyl phosphate       169.073            C(C(=O)COP(=O)(O)O)N
CID000000006    DNCB                                202.552            C1=CC(=C(C=C1[N+](=O)[O-])[N+](=O)[O-])Cl
CID000000007    9-ethyladenine                      163.18             CCN1C=NC2=C1N=CN=C2N
```

### <a name="OtherFiles"></a> 7. Other Files

#### 1. Other File 1

Text.

#### 2. Other File 2

Text

### <a name="Comments"></a> 8. Comments

1. <a name="noteCIDvsSID"></a>The prefix `CID` actually stands for "Compound ID", and should be distinguinshed from another chemical ID prefix (not in this dataset), `SID`, which stands for "Substance ID". These are PubChem conventions.
1. <a name="noteFileLegibility"></a> The file is not spaced as shown. The file is tab-delimited or space delimited. The spacing was added for better legibility.
