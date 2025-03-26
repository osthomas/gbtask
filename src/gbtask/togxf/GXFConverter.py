import logging
from collections import OrderedDict

import networkx as nx
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from gbtask.utils.TypeTree import DefaultTypeTree, TypeTree

log = logging.getLogger(__name__)


class GXFConverter:
    def __init__(self, record: SeqRecord, type_tree: TypeTree = DefaultTypeTree):
        self.record: SeqRecord = record
        self.type_tree: TypeTree = type_tree
        self._no_compound: set[str] = self._get_no_compound()
        self.format: str | None = None
        self.source: str = "gbtask_gxfconverter"

    def _get_no_compound(self):
        """
        Walk the underlying type tree to identify feature types that should not
        be splittable in the output format, ie. they should always cover their
        whole range, even if they are features with a CompoundLocation that
        have gaps.
        Anything at gtf_level 'transcript_id' or above in the hierarchy is not
        splittable.
        """
        # Find nodes at gtf_level transcript_id
        nodes = []
        for node, gtf_level in self.type_tree._tree.nodes(data="gtf_level"):
            if gtf_level == "transcript_id":
                nodes.append(node)
        # Find all their parent nodes
        unsplittable = set()
        for n in nodes:
            for x in nx.dfs_tree(self.type_tree._tree.reverse(), n).nodes:
                unsplittable.add(x)
        return unsplittable

    def convert_type(self, type: str):
        if not self.format:
            return type
        return self.type_tree.alias(type, self.format)

    def _feature_locations(self, feature: SeqFeature):
        """
        Return a list of one or more 3-tuples of (start, end, strand) for `feature`.
        For SimpleLocation features, the returned list will always be length 1.
        For CompoundLocation features, the length depends on whether the feature
        is splittable in the output or not:
        If it is, one tuple is returned per part.
        If it is not, the location returned covers the whole feature range without gaps.
        """
        if not feature.location:
            return [(".", ".", ".")]
        locs = []
        if feature.type not in self._no_compound:
            for part in feature.location.parts:
                start = part.start
                end = part.end  # type: ignore
                strand = part.strand
                locs.append((start, end, strand))
        else:
            locs.append(
                (
                    feature.location.start,
                    feature.location.end,
                    feature.location.strand,
                )
            )
        # Convert all locations to strings
        out = []
        for loc in locs:
            start, end, strand = loc
            # GFF locations are 1-based and end-inclusive
            start = int(start + 1)  # type: ignore
            start = str(start)
            end = str(int(end))  # type: ignore
            if strand == 1:
                strand = "+"
            elif strand == -1:
                strand = "-"
            elif strand == 0:
                strand = "?"
            else:
                strand = "."
            out.append((start, end, strand))
        return out

    def _feature_phases(self, feature: SeqFeature):
        """
        Return phases for each part of a feature.
        """
        if feature.type != "CDS":
            nparts = len(feature.location.parts) if feature.location else 1
            return ["."] * nparts
        cs = feature.qualifiers.get("codon_start", [1])[0]
        phase = int(cs) - 1  # codon_start is 1,2,3; phase is 0,1,2
        phases = [str(phase)]
        if not feature.location:
            return phases
        for part in feature.location.parts[:-1]:
            # for multi-part CDS: we have the first phase and need to calculate the
            # subsequent phases
            next = next_phase(int(phases[-1]), len(part))
            phases.append(str(next))
        return phases

    def feature_to_lines(self, feature: SeqFeature):
        raise NotImplementedError

    def header(self):
        raise NotImplementedError

    def convert(self):
        yield self.header()
        for feature in self.record.features:
            try:
                yield self.feature_to_lines(feature)
            except Exception as e:
                log.error(str(e))


class GFFConverter(GXFConverter):
    def __init__(self, record: SeqRecord, type_tree: TypeTree = DefaultTypeTree):
        super().__init__(record, type_tree)
        self.format: str | None = "gff"
        self.source: str = "gbtask_togff"

    def header(self):
        return f"##gff-version 3\n##sequence-region {self.record.name} 1 {len(self.record)}"

    def _feature_attributes(self, feature: SeqFeature):
        """
        Return qualifiers of `feature` formatted as GFF attributes.
        """
        attributes = OrderedDict()
        if "ID" not in feature.qualifiers or len(feature.qualifiers["ID"]) > 1:
            raise KeyError(
                "Feature does not have an ID qualifier. Run through geneify first!"
            )
        # ID and Parent first
        qid = feature.qualifiers["ID"][0]
        attributes["ID"] = qid
        qparents = ",".join(feature.qualifiers.get("Parent", []))
        if qparents:
            attributes["Parent"] = qparents
        # Insert all other qualifiers as attributes
        for k, v in feature.qualifiers.items():
            if k in ("ID", "Parent"):
                # already inserted
                continue
            attributes[k] = ",".join(v)
        # Merge individual attribute key-value pairs and concatenate all
        attributes = ";".join([f"{k}={v}" for k, v in attributes.items()])
        return attributes

    def feature_to_lines(self, feature: SeqFeature):
        chr = self.record.name
        source = self.source
        score = "."
        type = self.convert_type(feature.type)
        locs = self._feature_locations(feature)
        phases = self._feature_phases(feature)
        attributes = self._feature_attributes(feature)
        # All fields are identical across all parts of a compound feature,
        # except for location and phase
        lines = []
        for loc, phase in zip(locs, phases):
            start, end, strand = loc
            line = "\t".join(
                [chr, source, type, start, end, score, strand, phase, attributes]
            )
            lines.append(line)
        return "\n".join(lines)


class GTFConverter(GXFConverter):
    def __init__(self, record: SeqRecord, type_tree: TypeTree = DefaultTypeTree):
        super().__init__(record, type_tree)
        self.format: str | None = "gtf"
        self.source: str = "gbtask_togtf"

    def header(self):
        return "##format: gtf"

    def _feature_attributes(self, feature: SeqFeature):
        """
        Return qualifiers of `feature` formatted as GTF attributes.
        """
        attributes = OrderedDict()
        if (
            "gene_id" not in feature.qualifiers
            or len(feature.qualifiers["gene_id"]) > 1
        ):
            raise KeyError(
                "Feature does not have a single gene_id qualifier. A unique gene_id is required for GTF."
            )
        # gene_id and transcript_id first
        qgid = feature.qualifiers["gene_id"][0]
        attributes["gene_id"] = [qgid]
        gtxid = feature.qualifiers.get("transcript_id", [None])[0]
        if gtxid:
            attributes["transcript_id"] = [gtxid]
        for k, v in feature.qualifiers.items():
            # flatten all key-value pairs for one key.
            # Keys are repeated in GTF, it's not enough to just separate values
            # by comma
            attributes[k] = "; ".join([f'{k} "{v_}"' for v_ in v])
        # Merge all attributes
        attributes = "; ".join([f"{v}" for v in attributes.values()])
        return attributes

    def _txid_feature_lines(self, txid_feature: SeqFeature):
        """
        Generate GTF lines for one feature with one transcript_id, potentially
        expanding it over multiple lines if it is a compound feature.
        """
        lines = []
        chr = self.record.name
        source = self.source
        score = "."
        type = self.convert_type(txid_feature.type)
        locs = self._feature_locations(txid_feature)
        phases = self._feature_phases(txid_feature)
        attributes = self._feature_attributes(txid_feature)
        # All fields are identical across all parts of a compound feature,
        # except for location and phase
        for loc, phase in zip(locs, phases):
            start, end, strand = loc
            line = "\t".join(
                [chr, source, type, start, end, score, strand, phase, attributes]
            )
            lines.append(line)
        return "\n".join(lines)

    def feature_to_lines(self, feature: SeqFeature):
        """
        In GTF, a feature cannot be part of multiple 'transcript_id's.
        A feature that is part of more than one 'transcript_id' should be
        repeated on multiple lines.
        """
        txids = feature.qualifiers.get("transcript_id", [None])
        lines = []
        for txid in txids:
            # Copy feature, keeping all information the same, except for the
            # transript_id (if one is available).
            txid_feature = SeqFeature()
            txid_feature.id = feature.id
            txid_feature.type = feature.type
            txid_feature.location = feature.location
            for k, v in feature.qualifiers.items():
                if k == "transcript_id":
                    txid_feature.qualifiers[k] = [txid]
                else:
                    txid_feature.qualifiers[k] = v
            lines.append(self._txid_feature_lines(txid_feature))
        return "\n".join(lines)


def next_phase(phase: int, length: int):
    """
    Return phase of next CDS based on phase and length of current CDS.

    Examples
    --------
    [X] indicates first nucleotide of codon

    [A]TG---[G]AG
    >>> next_phase(0, 3)
    0

    A[T]G---G[A]G
    >>> next_phase(1, 3)
    1

    AT[G]---GA[G]
    >>> next_phase(2, 3)
    2

    [G]AT[G]---GA[G]
    >>> next_phase(0, 4)
    2

    G[G]AT---[G]AG
    >>> next_phase(1, 4)
    0

    GG[A]T---G[A]G
    >>> next_phase(2, 4)
    1

    G[G]AT[G]---GA[G]
    >>> next_phase(1, 5)
    2
    """
    # Make sure phase is in [0, 2]
    phase = phase % 3
    # nucleotides at end of first CDS that are part of an incomplete codon
    extra = (length - phase) % 3
    return (3 - extra) % 3
