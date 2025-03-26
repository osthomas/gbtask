from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbtask.utils import futils


def splice(record: SeqRecord, exon: str):
    # "Splicing" is done by concatening slices of the SeqRecord that correspond
    # to contiguous regions covered by exon features.
    # To properly handle all features, we need to do a bit more work:
    # Slicing only keeps features that are fully contained within the slice,
    # there is no automatic truncation.

    # First, make SimpleLocation features from CompoundLocation features
    record.features = _decompound(record.features)

    # Collect exons and determine contiguous intervals, to handle cases like:
    # exon1 exon2
    # vvvv  vvvv
    # ----  ----
    # ----------
    # ^^^^^^^^^^
    #   exon3
    exons = [f for f in record.features if f.type == exon]
    intervals = []
    for feature in record.features:
        if feature.type != exon or not feature.location:
            continue
        intervals.append((feature.location.start, feature.location.end))
    intervals = futils.reduce_intervals(intervals)
    if not intervals:
        # There are no exons with proper locations
        intervals = [(0, 0)]
    # create a fake spanning feature for each interval to use for truncation
    exon_spans = [SeqFeature(SimpleLocation(i[0], i[1])) for i in intervals]
    # Truncate features to the range of the exon slices. This ensures that they
    # will be retained when the SeqRecord is sliced.
    truncated_features = []
    exons = [f for f in record.features if f.type == exon]
    not_exons = [f for f in record.features if f.type != exon]
    for exon_span in exon_spans:
        for feature in not_exons:
            to_add = _truncate(feature, to=exon_span)
            if to_add:
                truncated_features.append(to_add)
    record.features = exons + truncated_features

    # Create a slice for each exon, then concatenate slices
    slices = [slice(i[0], i[1]) for i in intervals]
    segments = [record[slice] for slice in slices]
    spliced = segments.pop(0)
    while segments:
        spliced += segments.pop(0)
    # TODO ?  merge features which spanned multiple exons before and are now
    # split, but with adjacent ends.
    return spliced


def _decompound(features: list[SeqFeature]):
    """
    Slicing SeqRecords with SeqFeatures only keeps features which are fully
    contained in the slice.
    CompoundLocation features are considered whole, not per part. Therefore,
    split CompoundLocation features into SimpleLocation features.
    """
    out = []
    for feature in features:
        if not feature.location:
            continue
        for part in feature.location.parts:
            part_f = futils.copy_feature(feature)
            part_f.location = part
            out.append(part_f)
    return out


def _truncate(feature: SeqFeature, to: SeqFeature, stranded: bool = False):
    """
    Limit the range of `feature` to the range of feature `to`.


    Parameters
    ----------

    feature
        The feature to truncate
    to
        The feature defining the location limits beyond which truncation is applied
    stranded
        If True, only considers `features` that are on the same strand as `to`.

    Returns
    -------
    SeqFeature if `feature` can be truncated to `to`, None otherwise.
    """
    new = None
    if not to.location or not feature.location:
        return new
    if stranded and feature.location.strand != to.location.strand:
        return new
    if set(feature.location).intersection(to.location):
        # TODO: more efficient implementation of interval intersection
        new = futils.copy_feature(feature)
        if not isinstance(new.location, SimpleLocation):
            raise TypeError("Only SimpleLocation is supported")
        new_start = (
            to.location.start
            if new.location.start < to.location.start
            else new.location.start
        )
        new_end = (
            to.location.end if new.location.end > to.location.end else new.location.end
        )
        new.location = SimpleLocation(new_start, new_end, feature.location.strand)
    return new
