"""
Helpers to use in tests
"""

from Bio.SeqFeature import SeqFeature, SimpleLocation

# Wrappers to easily generate SeqFeatures.


def feature(
    type: str = "notype", start: int = 0, end: int = 10, strand=None, id="", **kwargs
):
    f = SeqFeature(
        SimpleLocation(start, end, strand), id=id, type=type, qualifiers=kwargs
    )
    return f


def mrna(start: int = 0, end: int = 10, strand=None, id="", **kwargs):
    return feature("mRNA", start, end, strand, id=id, **kwargs)


def exon(start: int = 0, end: int = 10, strand=None, id="", **kwargs):
    return feature("exon", start, end, strand, id=id, **kwargs)


def cds(start: int = 0, end: int = 10, strand=None, id="", **kwargs):
    return feature("CDS", start, end, strand, id=id, **kwargs)


def gene(name: str, start: int = 0, strand=None, id="", **kwargs):
    g = [
        feature("gene", start, start + 100, strand, id, **kwargs),
        mrna(start + 10, start + 90, strand, id, **kwargs),
        exon(start + 20, start + 30, strand, id, **kwargs),
        exon(start + 60, start + 80, strand, id, **kwargs),
        cds(start + 25, start + 30, strand, id, **kwargs),
    ]
    for f in g:
        f.qualifiers["gene"] = [name]
    return g


def new_gene(name, start):
    fdict = {
        "gene" + name: feature("gene", start + 0, start + 100, gene=["gene" + name]),
        "mrna" + name: feature("mRNA", start + 0, start + 100, gene=["gene" + name]),
        "exon1" + name: feature("exon", start + 10, start + 40, gene=["gene" + name]),
        "exon2" + name: feature("exon", start + 50, start + 60, gene=["gene" + name]),
        "cds" + name: feature("CDS", start + 51, start + 55, gene=["gene" + name]),
    }
    for key, f in fdict.items():
        f.id = key
    return fdict


def features():
    features = {}
    geneA = new_gene("A", 0)
    features.update(geneA)
    # Gene B, same location as Gene A
    geneB = new_gene("B", 0)
    features.update(geneB)
    # Copy of Gene A at another locus
    geneA2 = new_gene("A2", 300)
    for f in geneA2.values():
        f.qualifiers["gene"] = ["geneA"]
    features.update(geneA2)
    # Dangling features
    features.update({"dangling1": feature("misc_feature", start=0, end=100)})
    features.update({"dangling2": feature("regulatory", start=500, end=600)})
    for key, f in features.items():
        f.id = key
    return features
