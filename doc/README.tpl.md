# gbtask: CLI tool for manipulation of Genbank files

![Badge: pytest](https://github.com/osthomas/gbtask/actions/workflows/pytest.yaml/badge.svg)

$if(toc)$
$table-of-contents$

$endif$

`gbtask` is a set of commands for manipulation of Genbank files built on top
of Biopython.

Supported commands are:

* **geneify**: group features by qualifiers and overlap and build a feature hierarchy for export to GFF/GTF
* **modify**: padding, add/remove features, set qualifiers and annotation
* **splice**: carve out exon sequences and their overlapping features
* **togff/togtf**: leverage `geneify`'s hierarchization and export features to G(T/F)F
* **tofasta**: convert sequences to FASTA format

Its intended use case is to handle Genbank files storing reference information for specific genes of interest or plasmid annotations.
These potentially require occasional updates and are often not directly amenable to downstream tools which rely on GTF/GFF annotation.


## Installation

```
pip install git+https://www.github.com/osthomas/gbtask
```


## Examples

### Build a reference for a specific gene of interest

Start with an empty Genbank file:

```
${example_cat.sh()}


${example_cat.sh.out()}

```


Add exons:

```
${example_add1.sh()}


${example_add1.sh.out()}

```


Pad sequence, add exons and CDS, assign a custom grouping qualifier to all features:

```
${example_add2.sh()}


${example_add2.sh.out()}

```

Perform geneification:

```
${example_geneify.sh()}


${example_geneify.sh.out()}

```

Convert to GTF:

```
${example_togtf.sh()}


${example_togtf.sh.out()}

```

Convert to GFF:

```
${example_togff.sh()}


${example_togff.sh.out()}

```

Splice:

```
${example_splice.sh()}


${example_splice.sh.out()}

```

## geneify

`geneify` groups features by one or more qualifiers and by overlap, builds a hierarchy of features according to their types, and adds new qualifiers to each feature to bake the hierarchy into the Genbank format for subsequent export to GTF or GFF (qualifiers `ID` and `Parent` for GFF, `gene_id` and `transcript_id` for GTF).

`geneify` also infers exons (based on gapped mRNA features), mRNAs (based on exons) and UTRs (based on exon and CDS features).

This was inspired by BioPerl's
[Unflattener](https://metacpan.org/dist/BioPerl/source/lib/Bio/SeqFeature/Tools/Unflattener.pm)
and
[bp_genbank2gff3](https://metacpan.org/dist/BioPerl/view/bin/bp_genbank2gff3).

The feature hierarchy for supported features:

![Default Type Tree](doc/typetree.png)

Feature types that are not part of the hierarchy are assigned to the root feature type (gene).
Missing intermediates are skipped.

```
${geneify_help.sh()}


${geneify_help.sh.out()}

```


## modify

```
${modify_help.sh()}


${modify_help.sh.out()}

```


## splice

Note that features which are split by the splicing operation (eg. `gene`) are
duplicated and occur once for every exon by which they are covered.

```
${splice_help.sh()}


${splice_help.sh.out()}

```

## togtf


```
${togtf_help.sh()}


${togtf_help.sh.out()}

```



## togff


```
${togff_help.sh()}


${togff_help.sh.out()}

```
