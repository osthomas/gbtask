# gbtask: CLI tool for manipulation of Genbank files

- [gbtask: CLI tool for manipulation of Genbank
  files](#gbtask-cli-tool-for-manipulation-of-genbank-files)
  - [Installation](#installation)
  - [Examples](#examples)
    - [Build a reference for a specific gene of
      interest](#build-a-reference-for-a-specific-gene-of-interest)
  - [geneify](#geneify)
  - [modify](#modify)
  - [togtf](#togtf)
  - [togff](#togff)

`gbtask` is a set of commands for manipulation of Genbank files built on top
of Biopython.

Supported commands are:

* **geneify**: group features by qualifiers and overlap and build a feature hierarchy for export to GFF/GTF
* **modify**: padding, add/remove features, set qualifiers and annotation
* **togff/togtf**: leverage `geneify`'s hierarchization and export features to G(T/F)F

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
cat example.gb

LOCUS       example                 100 bp     DNA     linear   UNA 01-JAN-2000
DEFINITION  natural linear DNA.
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      natural DNA sequence
  ORGANISM  unspecified
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 atgcatgcat atgcatgcat atgcatgcat atgcatgcat atgcatgcat atgcatgcat
       61 atgcatgcat atgcatgcat atgcatgcat atgcatgcat
//
```


Add exons:

```
gbtask modify --feature 0 10 1 exon --feature -0 -10 1 exon -i example.gb > ex1.gb
cat ex1.gb

LOCUS       example                  100 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  natural linear DNA.
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      natural DNA sequence
  ORGANISM  unspecified
            .
FEATURES             Location/Qualifiers
     exon            1..10
     exon            91..100
ORIGIN
        1 atgcatgcat atgcatgcat atgcatgcat atgcatgcat atgcatgcat atgcatgcat
       61 atgcatgcat atgcatgcat atgcatgcat atgcatgcat
//
```


Pad sequence, add exons and CDS, assign a custom grouping qualifier to all features:

```
gbtask modify --pad-left 100 --pad-right AAAAGGGGGCCCCCTTTTT \
    --feature 10 20 1 exon --feature 11 19 1 CDS \
    -q all mygroupingqualifier mygene \
    -i ex1.gb > ex2.gb
cat ex2.gb

LOCUS       example                  219 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  natural linear DNA.
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      natural DNA sequence
  ORGANISM  unspecified
            .
FEATURES             Location/Qualifiers
     exon            101..110
                     /mygroupingqualifier="mygene"
     exon            191..200
                     /mygroupingqualifier="mygene"
     exon            11..20
                     /mygroupingqualifier="mygene"
     CDS             12..19
                     /mygroupingqualifier="mygene"
ORIGIN
        1 nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn
       61 nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn atgcatgcat atgcatgcat
      121 atgcatgcat atgcatgcat atgcatgcat atgcatgcat atgcatgcat atgcatgcat
      181 atgcatgcat atgcatgcat aaaagggggc ccccttttt
//
```

Perform geneification:

```
gbtask geneify --groupby mygroupingqualifier -i ex2.gb > geneified.gb
cat geneified.gb

LOCUS       example                  219 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  natural linear DNA.
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      natural DNA sequence
  ORGANISM  unspecified
            .
FEATURES             Location/Qualifiers
     exon            101..110
                     /mygroupingqualifier="mygene"
                     /ID="gene.01:exon.02"
                     /Parent="gene.01:mRNA.01"
                     /transcript_id="gene.01:mRNA.01"
                     /gene_id="gene.01"
     exon            191..200
                     /mygroupingqualifier="mygene"
                     /ID="gene.01:exon.03"
                     /Parent="gene.01:mRNA.01"
                     /transcript_id="gene.01:mRNA.01"
                     /gene_id="gene.01"
     exon            11..20
                     /mygroupingqualifier="mygene"
                     /ID="gene.01:exon.01"
                     /Parent="gene.01:mRNA.01"
                     /transcript_id="gene.01:mRNA.01"
                     /gene_id="gene.01"
     CDS             12..19
                     /mygroupingqualifier="mygene"
                     /ID="gene.01:CDS.01"
                     /Parent="gene.01:mRNA.01"
                     /transcript_id="gene.01:mRNA.01"
                     /gene_id="gene.01"
     mRNA            join(11..20,101..110,191..200)
                     /ID="gene.01:mRNA.01"
                     /Parent="gene.01"
                     /transcript_id="gene.01:mRNA.01"
                     /gene_id="gene.01"
     5'UTR           11
                     /ID="gene.01:5'UTR.01"
                     /Parent="gene.01:mRNA.01"
                     /transcript_id="gene.01:mRNA.01"
                     /gene_id="gene.01"
     3'UTR           join(20,101..110,191..200)
                     /ID="gene.01:3'UTR.01"
                     /Parent="gene.01:mRNA.01"
                     /transcript_id="gene.01:mRNA.01"
                     /gene_id="gene.01"
     gene            11..200
                     /ID="gene.01"
                     /gene_id="gene.01"
ORIGIN
        1 nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn
       61 nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn nnnnnnnnnn atgcatgcat atgcatgcat
      121 atgcatgcat atgcatgcat atgcatgcat atgcatgcat atgcatgcat atgcatgcat
      181 atgcatgcat atgcatgcat aaaagggggc ccccttttt
//
```

Convert to GTF:

```
gbtask togtf -i geneified.gb

##format: gtf
example	gbtask_togtf	exon	101	110	.	+	.	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; mygroupingqualifier "mygene"; ID "gene.01:exon.02"; Parent "gene.01:mRNA.01"
example	gbtask_togtf	exon	191	200	.	+	.	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; mygroupingqualifier "mygene"; ID "gene.01:exon.03"; Parent "gene.01:mRNA.01"
example	gbtask_togtf	exon	11	20	.	+	.	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; mygroupingqualifier "mygene"; ID "gene.01:exon.01"; Parent "gene.01:mRNA.01"
example	gbtask_togtf	CDS	12	19	.	+	0	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; mygroupingqualifier "mygene"; ID "gene.01:CDS.01"; Parent "gene.01:mRNA.01"
example	gbtask_togtf	mRNA	11	200	.	+	.	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; ID "gene.01:mRNA.01"; Parent "gene.01"
example	gbtask_togtf	UTR	11	11	.	+	.	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; ID "gene.01:5'UTR.01"; Parent "gene.01:mRNA.01"
example	gbtask_togtf	UTR	20	20	.	+	.	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; ID "gene.01:3'UTR.01"; Parent "gene.01:mRNA.01"
example	gbtask_togtf	UTR	101	110	.	+	.	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; ID "gene.01:3'UTR.01"; Parent "gene.01:mRNA.01"
example	gbtask_togtf	UTR	191	200	.	+	.	gene_id "gene.01"; transcript_id "gene.01:mRNA.01"; ID "gene.01:3'UTR.01"; Parent "gene.01:mRNA.01"
example	gbtask_togtf	gene	11	200	.	+	.	gene_id "gene.01"; ID "gene.01"
```

Convert to GFF:

```
gbtask togff -i geneified.gb

##gff-version 3
##sequence-region example 1 219
example	gbtask_togff	exon	101	110	.	+	.	ID=gene.01:exon.02;Parent=gene.01:mRNA.01;mygroupingqualifier=mygene;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	exon	191	200	.	+	.	ID=gene.01:exon.03;Parent=gene.01:mRNA.01;mygroupingqualifier=mygene;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	exon	11	20	.	+	.	ID=gene.01:exon.01;Parent=gene.01:mRNA.01;mygroupingqualifier=mygene;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	CDS	12	19	.	+	0	ID=gene.01:CDS.01;Parent=gene.01:mRNA.01;mygroupingqualifier=mygene;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	mRNA	11	200	.	+	.	ID=gene.01:mRNA.01;Parent=gene.01;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	five_prime_utr	11	11	.	+	.	ID=gene.01:5'UTR.01;Parent=gene.01:mRNA.01;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	three_prime_utr	20	20	.	+	.	ID=gene.01:3'UTR.01;Parent=gene.01:mRNA.01;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	three_prime_utr	101	110	.	+	.	ID=gene.01:3'UTR.01;Parent=gene.01:mRNA.01;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	three_prime_utr	191	200	.	+	.	ID=gene.01:3'UTR.01;Parent=gene.01:mRNA.01;transcript_id=gene.01:mRNA.01;gene_id=gene.01
example	gbtask_togff	gene	11	200	.	+	.	ID=gene.01;gene_id=gene.01
```


## geneify

`geneify` groups features by one or more qualifiers and by overlap, builds a hierarchy of features according to their types, and adds new qualifiers to each feature to bake the hierarchy into the Genbank format for subsequent export to GTF or GFF (qualifiers `ID` and `Parent` for GFF, `gene_id` and `transcript_id` for GTF).

`geneify` also infers exons (based on gapped mRNA features), mRNAs (based on exons) and UTRs (based on exon and CDS features).

This was inspired by BioPerl's
[Unflattener](https://metacpan.org/dist/BioPerl/source/lib/Bio/SeqFeature/Tools/Unflattener.pm)
and
[bp_genbank2gff3](https://metacpan.org/dist/BioPerl/view/bin/bp_genbank2gff3).

The feature hierarchy for supported features:

![Default Type Tree](typetree.png)

Feature types that are not part of the hierarchy are assigned to the root feature type (gene).
Missing intermediates are skipped.

```
gbtask geneify -h

usage: gbtask geneify [-h] [-i INFILE] [-o OUTFILE] [--groupby [GROUPBY ...]]
                      [--groupby-overlap | --no-groupby-overlap]
                      [--infer-exons | --no-infer-exons]
                      [--infer-utrs | --no-infer-utrs]

options:
  -h, --help            show this help message and exit
  -i, --infile INFILE   Input file (default: <_io.TextIOWrapper name='<stdin>'
                        mode='r' encoding='utf-8'>)
  -o, --outfile OUTFILE
                        Output file (default: <_io.TextIOWrapper
                        name='<stdout>' mode='w' encoding='utf-8'>)
  --groupby [GROUPBY ...]
                        Names of qualifiers for first-level grouping of
                        features (default: ['gene'])
  --groupby-overlap, --no-groupby-overlap
                        Group features by overlap after grouping by
                        qualifiers? Warning, this may lead to undesirable
                        separation of features if there is no scaffolding
                        super-feature that covers all sub-features. (default:
                        False)
  --infer-exons, --no-infer-exons
                        Infer mRNA/exon features per feature group, depending
                        on what is available (default: True)
  --infer-utrs, --no-infer-utrs
                        Infer UTR features per feature group based on exons
                        and CDS (default: True)
```


## modify

```
gbtask modify -h

usage: gbtask modify [-h] [-i INFILE] [-o OUTFILE] [-a key value]
                     [-q type key value] [-f start end strand type]
                     [--include [INCLUDE ...] | --exclude [EXCLUDE ...]]
                     [--pad-left PAD_LEFT] [--pad-right PAD_RIGHT]

options:
  -h, --help            show this help message and exit
  -i, --infile INFILE   Input file (default: <_io.TextIOWrapper name='<stdin>'
                        mode='r' encoding='utf-8'>)
  -o, --outfile OUTFILE
                        Output file (default: <_io.TextIOWrapper
                        name='<stdout>' mode='w' encoding='utf-8'>)
  -a, --annotation key value
                        Set annotation 'key' to 'value' for all records. Can
                        be given multiple times. Example: -a organism 'mus
                        musculus' (default: [])
  -q, --qualifier type key value
                        For all features of type 'type', set qualifier 'key'
                        to 'value'. Can be given multiple times. The special
                        value 'all' for 'type' can be given to affect all
                        features.Example: -q all locus_tag geneA (default: [])
  -f, --feature start end strand type
                        Add a new feature. Can be given multiple times. NOTE:
                        Locations 'start' and 'end' are 0-based and end-
                        exclusive! If they start with '-', counting begins
                        from the end of the sequence. '0' refers to the first,
                        '-0' refers to the last base of the record. (default:
                        [])

Include/Exclude Features:
  Feature types to include or exclude. Only --include or --exclude can be
  specified, not both. By default, all feature types are included.

  --include [INCLUDE ...]
                        Only output these feature types. (default: [])
  --exclude [EXCLUDE ...]
                        Exclude these feature types from the output. (default:
                        [])

Padding:
  --pad-left PAD_LEFT   Pad sequence to the left. This can either be a number
                        to add this many 'N' characters, or a string to pad
                        with a specific sequence. (default: None)
  --pad-right PAD_RIGHT
                        Like --pad-left, but for the right side of the
                        sequence. (default: None)
```


## togtf


```
gbtask togtf -h

usage: gbtask togtf [-h] [-i INFILE] [-o OUTFILE]

options:
  -h, --help            show this help message and exit
  -i, --infile INFILE   Input file (default: <_io.TextIOWrapper name='<stdin>'
                        mode='r' encoding='utf-8'>)
  -o, --outfile OUTFILE
                        Output file (default: <_io.TextIOWrapper
                        name='<stdout>' mode='w' encoding='utf-8'>)
```



## togff


```
gbtask togff -h

usage: gbtask togff [-h] [-i INFILE] [-o OUTFILE]

options:
  -h, --help            show this help message and exit
  -i, --infile INFILE   Input file (default: <_io.TextIOWrapper name='<stdin>'
                        mode='r' encoding='utf-8'>)
  -o, --outfile OUTFILE
                        Output file (default: <_io.TextIOWrapper
                        name='<stdout>' mode='w' encoding='utf-8'>)
```
