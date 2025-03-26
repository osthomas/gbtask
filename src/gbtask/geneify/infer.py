import logging
from typing import Literal

from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from gbtask.utils.futils import (
    common_strand,
    covered_intervals,
    logstr,
    reduce_intervals,
    strands_equal,
)

log = logging.getLogger(__name__)


class Inferrer:
    def __init__(self, features: list[SeqFeature]):
        raise NotImplementedError

    def infer(self) -> list[SeqFeature]:
        raise NotImplementedError


class ExonInferrer(Inferrer):
    def __init__(
        self, features: list[SeqFeature], exon: str = "exon", rna: str = "mRNA"
    ):
        """
        Infer exons from a list of features.

        How inference is performed depends on the feature types available in the
        feature list:

        1. If features of type `rna` are available, create new exons over all
           positions covered by the RNA, skipping gaps. Contiguous regions
           that are fully covered by existing features of type `exon` with the
           correct strand are skipped, ie. exons are not duplicated. The
           generated exons are added to the list of known exons, together with
           the provided exons.
        2. If features of type `exon` are available, but not of type `rna`, new
           RNA features are created by collecting exons per strand.

        Parameters
        ----------

        exon
            Name of the feature type for exons.
        rna
            Name of the feature type for RNAs, eg. "mRNA", "tRNA", ...
        """
        if not strands_equal(features):
            msg = "Strands of all features must be equal for exon inference. Split features by strand first."
            raise ValueError(msg)
        self.features: list[SeqFeature] = features
        self.exon: str = exon
        self.rna: str = rna
        self._new_exons: list[SeqFeature] = []  # storage space for inferred exons
        self._new_rnas: list[SeqFeature] = []  # storage space for inferred rnas

    @property
    def exons(self):
        """
        Return a list of all known exons, pre-existing and already inferred.
        """
        out = []
        for feature in self.features:
            if feature.type == self.exon:
                out.append(feature)
        out.extend(self._new_exons)
        return out

    @property
    def rnas(self):
        """
        Return a list of all known RNAs, pre-existing and already inferred.
        """
        out = []
        for feature in self.features:
            if feature.type == self.rna:
                out.append(feature)
        out.extend(self._new_rnas)
        return out

    def is_compatible(self, exon: SeqFeature, rna: SeqFeature):
        """
        Return True if the exon is compatible with the rna feature.
        Return False for incompatible features.
        """
        exon_loc = exon.location
        rna_loc = rna.location
        compat = True
        if exon_loc is None or rna_loc is None:
            compat = False
        elif exon_loc.strand != rna_loc.strand:
            compat = False
        elif exon_loc.start > rna_loc.end:
            compat = False
        elif exon_loc.end < rna_loc.start:
            compat = False
        compat = compat and set(rna).issuperset(set(exon))
        if not compat:
            log.warning(f"Exon {logstr(exon)} is incompatible with RNA {logstr(rna)}")
        return compat

    def _exons_from_rna(self, rna: SeqFeature):
        if not rna.location:
            return []
        exons = []
        covered = []
        # Collect positions that are already covered by exons
        for feature in self.exons:
            if self.is_compatible(feature, rna):
                feature_loc = feature.location
                if feature_loc:
                    covered.append((feature_loc.start, feature_loc.end))
        # Create new exons for stretches of RNA that are not covered by
        # existing exons
        for part in rna.location.parts:
            cover = (part.start, part.end)
            if cover not in covered:
                exon = SeqFeature(part, type="exon")
                exons.append(exon)
        return exons

    def _rna_from_exons(self, exons: list[SeqFeature]):
        rnas = []
        locations = covered_intervals(exons)
        if len(locations) == 1:
            loc = locations[0]
        else:
            loc = CompoundLocation(locations)
        rnas.append(SeqFeature(loc, type=self.rna))
        return rnas

    def infer(self):
        """
        Infer RNA and Exon features and return a list of new features.
        Pre-existing exons/RNAs are not returned.
        """
        if self.rnas:
            for rna in self.rnas:
                self._new_exons.extend(self._exons_from_rna(rna))
        elif self.exons:
            self._new_rnas.extend(self._rna_from_exons(self.exons))
        return self._new_exons + self._new_rnas


class UTRInferrer(Inferrer):
    def __init__(
        self,
        features: list[SeqFeature],
        exon: str = "exon",
        cds: str = "CDS",
        utr3p: str = "3'UTR",
        utr5p: str = "5'UTR",
        utr: str = "UTR",
    ):
        """
        Parameters
        ----------

        features
            List of features for which to infer UTR locations
        exon
            Name of the feature type for exons.
        cds
            Name of the feature type for CDS.
        utr3p
            Name of the feature type for 3' UTRs.
        utr5p
            Name of the feature type for 5' UTRs.
        utr
            Name of the feature type for an unspecified UTR (either 5' or 3')
        """
        if not strands_equal(features):
            msg = "Strands of all features must be equal for UTR inference. Split features by strand first."
            raise ValueError(msg)
        self.features: list[SeqFeature] = features
        self.exon: str = exon
        self.cds: str = cds
        self.utr3p: str = utr3p
        self.utr5p: str = utr5p
        self.utr: str = utr
        self._strand = common_strand(self.features)

    def _utr_intervals(self):
        # Collect relevant features with known locations
        left_utr_locs = []
        right_utr_locs = []
        exons = []
        cdss = []
        for feature in self.features:
            if feature.type == self.exon:
                exons.append(feature)
            elif feature.type == self.cds:
                cdss.append(feature)

        if not exons or not cdss:
            return (left_utr_locs, right_utr_locs)

        cds_range = [
            min([cds.location.start for cds in cdss]),
            max([cds.location.end for cds in cdss]),
        ]

        # Move through exons and collect contiguous intervals
        for exon in exons:
            cur_pos = exon.location.start  # snap to start of current exon
            if cur_pos < cds_range[0]:
                # in left UTR
                next_pos = min(cds_range[0], exon.location.end)
                if next_pos > cur_pos:
                    left_utr_locs.append((cur_pos, next_pos))
                cur_pos = next_pos
            if cur_pos >= cds_range[0]:
                # in CDS or in right UTR
                cur_pos = max(cur_pos, cds_range[1])  # snap to end of CDS
                if cur_pos > exon.location.end:
                    # Exon fully covers CDS, no right UTR
                    continue
                next_pos = exon.location.end
                if next_pos > cur_pos:
                    right_utr_locs.append((cur_pos, next_pos))
                cur_pos = next_pos

        # Reduce
        left_utr_locs = reduce_intervals(left_utr_locs)
        right_utr_locs = reduce_intervals(right_utr_locs)
        return left_utr_locs, right_utr_locs

    def _utr_feature(
        self, intervals: list[tuple[int, int]], side: Literal["left", "right"]
    ):
        locs = [SimpleLocation(i[0], i[1], self._strand) for i in intervals]
        if len(locs) > 1:
            loc = CompoundLocation(locs)
        else:
            loc = locs[0]
        utr_types = {"left": self.utr5p, "right": self.utr3p}
        if self._strand == 1:
            utr_type = utr_types[side]
        elif self._strand == -1:
            # swap sides for UTR on reverse strand
            if side == "left":
                side = "right"
            elif side == "right":
                side = "left"
            utr_type = utr_types[side]
        else:
            utr_type = self.utr
            log.warning("Strand of inferred UTRs is unclear, cannot determine 3' or 5'")
        utr = SeqFeature(loc, type=utr_type)
        return utr

    def _preexisting(self, new_utr: SeqFeature):
        """
        Return True if an inferred UTR already exists in the feature list.
        """
        for feature in self.features:
            same = True
            same = same and feature.type == new_utr.type
            same = same and feature.location == new_utr.location
            if same:
                return True
        return False

    def infer(self):
        left_locs, right_locs = self._utr_intervals()
        utrs = []
        if left_locs:
            left_utr = self._utr_feature(left_locs, "left")
            if not self._preexisting(left_utr):
                utrs.append(left_utr)
        if right_locs:
            right_utr = self._utr_feature(right_locs, "right")
            if not self._preexisting(right_utr):
                utrs.append(right_utr)
        return utrs
