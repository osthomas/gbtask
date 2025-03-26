"""
futils - feature utilities
Utility functions for handling SeqFeature objects.
"""

import logging
from typing import Iterable, Literal

from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

log = logging.getLogger(__name__)


def logstr(feature: SeqFeature, qualifiers: Iterable | None | Literal[True] = None):
    """
    Compact representation of feature for logging purposes.

    Parameters
    ----------

    qualifiers
        Which qualifiers to output.
        If this is a list, output qualifiers with these keys.
        If None, output no qualifiers.
        If True, output all qualifiers.
    """
    if qualifiers is None:
        qualifiers = []
    elif qualifiers is True:
        qualifiers = feature.qualifiers.keys()
    qvals = []
    for key in qualifiers:
        v = feature.qualifiers.get(key)
        v = str(v)
        if len(v) > 30:
            v = v[:30] + "..."
        qvals.append(f"{key}={v}")
    if qvals:
        qvals = ", ".join(qvals)
        return f"{feature.type} (ID: {feature.id}, location: {feature.location}, qualifiers: {qvals})"
    else:
        return f"{feature.type} (ID: {feature.id}, location: {feature.location})"


def combine_feature_locations(features: list[SeqFeature]):
    locs = []
    for feature in features:
        if not feature.location:
            continue
        locs.extend(feature.location.parts)
    if len(locs) > 1:
        return CompoundLocation(locs)
    else:
        return locs[0]


def covers(x: SeqFeature, y: SeqFeature, stranded: bool = False):
    """
    Return True if feature `x` fully covers all locations of feature `y`. Gaps
    in feature `x` do not contribute to coverage.

    Examples
    --------

                 1         2         3         4
        1        0         0         0         0
    f1  --------->
    f2                     ---------->
    f3                          <-----
    f4            -------->           -------->

    >>> from Bio.SeqFeature import SimpleLocation, CompoundLocation
    >>> f1 = SeqFeature(SimpleLocation(1, 10, 1))
    >>> f2 = SeqFeature(SimpleLocation(20, 30, 1))
    >>> f3 = SeqFeature(SimpleLocation(25, 30, -1))
    >>> f4 = SeqFeature(CompoundLocation([SimpleLocation(11, 19, 1), SimpleLocation(31, 40, 1)]))

    >>> covers(f1, f2)
    False
    >>> covers(f2, f3)
    True
    >>> covers(f2, f3, stranded=True)
    False
    >>> covers(f3, f2)
    False
    >>> covers(f4, f2)
    False
    """
    if not x.location or not y.location:
        # no information
        return False
    ok = True
    if stranded:
        ok = ok and x.location.strand == y.location.strand
    for loc in y.location:
        ok = ok and loc in x
    return ok


def common_strand(features: list[SeqFeature]) -> Literal[1, -1, None]:
    """
    Return the strand of all `features` if it is identical for all. This is
    `None` if the strand of all features is `None`!
    `None` is also returned if strands differ between features.

    Examples
    --------
    >>> minus = SeqFeature(SimpleLocation(0, 100, -1))
    >>> plus = SeqFeature(SimpleLocation(0, 100, 1))
    >>> unknown = SeqFeature(SimpleLocation(0, 100, None))

    >>> common_strand([plus, plus])
    1

    >>> common_strand([minus, minus])
    -1

    >>> common_strand([unknown, unknown]) is None
    True

    >>> common_strand([minus, plus]) is None
    True
    """
    loc = combine_feature_locations(features)
    return loc.strand


def strands_equal(features: list[SeqFeature]) -> bool:
    """
    Return True if the strand of all `features` are identical.

    Examples
    --------
    >>> minus = SeqFeature(SimpleLocation(0, 100, -1))
    >>> plus = SeqFeature(SimpleLocation(0, 100, 1))
    >>> unknown = SeqFeature(SimpleLocation(0, 100, None))

    >>> strands_equal([minus, minus])
    True

    >>> strands_equal([plus, minus])
    False

    >>> strands_equal([unknown, unknown])
    True
    """
    if not features:
        return True
    else:
        strand = features[0].location.strand if features[0].location else None
    for feature in features:
        fstrand = feature.location.strand if feature.location else None
        if strand != fstrand:
            return False
        else:
            strand = fstrand
    return True


def has_strand(features: list[SeqFeature]) -> bool:
    """
    Return True if all features have a strand definition that is not 'None'.

    Examples
    --------
    >>> minus = SeqFeature(SimpleLocation(0, 100, -1))
    >>> plus = SeqFeature(SimpleLocation(0, 100, 1))
    >>> unknown = SeqFeature(SimpleLocation(0, 100, None))

    >>> has_strand([])
    False

    >>> has_strand([plus, plus])
    True

    >>> has_strand([plus, minus])
    True

    >>> has_strand([plus, unknown])
    False

    >>> has_strand([minus, unknown])
    False
    """
    if not features:
        return False
    try:
        return all(
            [f.location.strand in (1, -1) for f in features]  # pyright: ignore[reportOptionalMemberAccess]
        )
    except AttributeError:
        # no location or strand defined
        return False


def create_spanning(type: str, features: list[SeqFeature]) -> SeqFeature:
    """
    Create a feature spanning the extent of other features. Gaps in `features`
    are fully covered, the resulting feature has no gaps. The strand of the
    output feature is inferred via :func:`~common_strand`.

    Params
    ------

    type
        feature type of output feature
    features
        list of features that should be spanned
    strand
        strand of the output feature


    Examples
    --------

    >>> from Bio.SeqFeature import SimpleLocation, CompoundLocation
    >>> l1 = SimpleLocation(0, 100)
    >>> l2 = SimpleLocation(200, 300)
    >>> f1 = SeqFeature(l1)

    >>> create_spanning("gene", [f1])
    SeqFeature(SimpleLocation(ExactPosition(0), ExactPosition(100)), type='gene')

    >>> f2 = SeqFeature(l2)
    >>> create_spanning("gene", [f1, f2])
    SeqFeature(SimpleLocation(ExactPosition(0), ExactPosition(300)), type='gene')

    >>> f3 = SeqFeature(CompoundLocation([l1, l2]))
    >>> create_spanning("gene", [f3])
    SeqFeature(SimpleLocation(ExactPosition(0), ExactPosition(300)), type='gene')
    """
    loc = combine_feature_locations(features)
    return SeqFeature(SimpleLocation(min(loc), max(loc) + 1, loc.strand), type=type)


def covered_intervals(features: list[SeqFeature]) -> list[SimpleLocation]:
    """
    Return a list of `SimpleLocation`, each describing a unique range that is
    covered (without gaps) by the `features`.

    >>> from Bio.SeqFeature import SeqFeature
    >>> from Bio.SeqFeature import SimpleLocation as SL
    >>> from Bio.SeqFeature import CompoundLocation as CL

    >>> rna = SeqFeature(CL([SL(0, 10), SL(20, 30), SL(30, 40)]))
    >>> covered_intervals([rna])
    [SimpleLocation(ExactPosition(0), ExactPosition(10)), SimpleLocation(ExactPosition(20), ExactPosition(40))]

    >>> exon1 = SeqFeature(SL(0, 10))
    >>> exon2 = SeqFeature(SL(20, 30))
    >>> covered_intervals([exon1, exon2])
    [SimpleLocation(ExactPosition(0), ExactPosition(10)), SimpleLocation(ExactPosition(20), ExactPosition(30))]
    """
    covered = []
    for feature in features:
        if not feature.location:
            continue
        for part in feature.location.parts:
            covered.append((part.start, part.end))
    covered = reduce_intervals(covered)
    return [SimpleLocation(s, e, common_strand(features)) for s, e in covered]


def reduce_intervals(intervals: list[tuple[int, int]]):
    """
    Reduce adjacent contiguous *right-open* intervals to one bigger interval.

    Return
    ------

    A list of reduced intervals.


    Examples
    --------

    >>> reduce_intervals([(0, 1), (1, 2), (2, 3)])
    [(0, 3)]

    Intervals are sorted by start position:

    >>> reduce_intervals([(0, 1), (2, 3), (1, 2)])
    [(0, 3)]

    >>> reduce_intervals([(0, 1), (1, 2), (3, 4)])
    [(0, 2), (3, 4)]

    >>> reduce_intervals([(0, 1)])
    [(0, 1)]

    >>> reduce_intervals([(0, 1), (2, 3)])
    [(0, 1), (2, 3)]

    >>> reduce_intervals([(0, 0), (1, 1)])
    [(0, 0), (1, 1)]

    Overlapping intervals are merged:

    >>> reduce_intervals([(0, 3), (2, 5)])
    [(0, 5)]

    >>> reduce_intervals([])
    []
    """
    if not intervals:
        return []
    # Sort by start position
    intervals = sorted(intervals, key=lambda x: x[0])
    new = []
    # Fetch one interval
    i1 = intervals.pop(0)
    while intervals:
        i2 = intervals.pop(0)
        fused = fuse_intervals(i1, i2)
        if len(fused) == 1:
            #  Intervals were fused, continue to checking next interval in list
            i1 = fused[0]
        elif len(fused) == 2:
            # Intervals were not fused, i1 is done.
            # Continue checking if i2 can be fused with next intervals
            new.append(fused[0])
            i1 = fused[1]
    new.append(i1)
    return new


def fuse_intervals(i1: tuple[int, int], i2: tuple[int, int]):
    """
    Fuse adjacent right-open intervals.

    Return
    ------

    A list of tuples. If `i1` and `i2` were fused, the list has length 1, else
    it has length 2.


    Examples
    --------
    >>> fuse_intervals((0, 1), (1, 2))
    [(0, 2)]

    >>> fuse_intervals((1, 1), (1, 2))
    [(1, 2)]

    >>> fuse_intervals((1, 0), (1, 2))
    Traceback (most recent call last):
    ...
    ValueError: Interval ends must be >= interval starts

    >>> fuse_intervals((0, 1), (2, 3))
    [(0, 1), (2, 3)]
    """
    if i1[0] > i1[1] or i2[0] > i2[1]:
        raise ValueError("Interval ends must be >= interval starts")
    # Sort by start position
    if i1[0] > i2[0]:
        i1, i2 = i2, i1
    if i1[1] >= i2[0]:
        new = [(i1[0], i2[1])]
    else:
        new = [i1, i2]
    return new


def split_by_strand(features: list[SeqFeature]):
    """
    Split a list of features into a dictionary, separating features by strand.
    Features with strand `None` are put under key `"None"` (as a string).

    Examples
    --------

    >>> minus = SeqFeature(SimpleLocation(0, 100, -1))
    >>> plus = SeqFeature(SimpleLocation(0, 100, 1))
    >>> unknown = SeqFeature(SimpleLocation(0, 100, None))

    >>> split_by_strand([minus, plus, unknown])
    {-1: [SeqFeature(SimpleLocation(ExactPosition(0), ExactPosition(100), strand=-1))], 1: [SeqFeature(SimpleLocation(ExactPosition(0), ExactPosition(100), strand=1))], 'None': [SeqFeature(SimpleLocation(ExactPosition(0), ExactPosition(100)))]}
    """
    d = {}
    for feature in features:
        strand = feature.location.strand if feature.location else None
        strand = strand or "None"
        if strand not in d:
            d[strand] = []
        d[strand].append(feature)
    return d


def copy_feature(feature: SeqFeature):
    """
    Return a copy of a SeqFeature.
    The qualifier dict of the new feature is a shallow copy of the original.
    """
    new = SeqFeature()
    new.location = feature.location
    new.type = feature.type
    new.id = feature.id
    # shallow copy of qualifier dict
    new.qualifiers = {k: v for k, v in feature.qualifiers.items()}
    return new
