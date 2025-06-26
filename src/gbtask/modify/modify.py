from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbtask.utils.futils import reduce_intervals


def modify(
    record: SeqRecord,
    pad_left: str,
    pad_right: str,
    include: list[str],
    exclude: list[str],
    features: list[tuple[str, str, str, str]],
    annotations: list[tuple[str, str]],
    qualifiers: list[tuple[str, str, str]],
    expansions: list[tuple[str, str]],
):
    if pad_left:
        record = pad2seq(pad_left) + record
    if pad_right:
        record = record + pad2seq(pad_right)
    record.features = [f for f in record.features if included(f, include, exclude)]  # type: ignore
    if features:
        add_features(record, features)
    if annotations:
        set_annotations(record, annotations)
    if qualifiers:
        set_qualifiers(record, qualifiers)
    if expansions:
        expand_features(record, expansions)
    return record


def pad2seq(pad: str):
    try:
        seq = Seq(int(pad) * "N")
    except:
        seq = Seq(pad)
    return seq


def included(feature: SeqFeature, include: list[str], exclude: list[str]):
    """
    Return True if a feature is included for processing.

    Inclusion takes precedence over exclusion is both `include` and `exclude`
    are specified.

    Examples
    --------
    >>> included(SeqFeature(type = "exon"), include = ["exon"], exclude = [])
    True
    >>> included(SeqFeature(type = "gene"), include = ["exon"], exclude = [])
    False
    >>> included(SeqFeature(type = "exon"), include = [], exclude = ["exon"])
    False
    >>> included(SeqFeature(type = "exon"), include = ["exon"], exclude = ["exon"])
    True

    Types are not case sensitive:
    >>> included(SeqFeature(type = "Cds"), include = ["cdS"], exclude = [])
    True
    """
    # do not check type names with case sensitivity
    include = [i.lower() for i in include]
    exclude = [e.lower() for e in exclude]
    ftype = feature.type.lower()
    if not include and not exclude:
        # include everything
        return True
    if include and ftype in include:
        return True
    if exclude and ftype not in exclude:
        return True
    return False


def set_annotations(record: SeqRecord, annotation_spec: list[tuple[str, str]]):
    """
    Set annotations for a SeqRecord in place.

    Parameters
    ----------
    annotation_spec
        A list of tuples.
        Each tuple `t` has 2 string elements.
        `t[0]` gives the key of the annotation.
        `t[1]` gives the values the annotation is set to.
    """
    for k, v in annotation_spec:
        # SeqRecord annotation keys are lowercase
        k = k.lower()
        if k == "locus":
            record.name = v
        else:
            record.annotations[k] = v
    return record


def set_qualifiers(record: SeqRecord, qualifier_spec: list[tuple[str, str, str]]):
    """
    Modify qualifier values for features in `record` in place.

    Parameters
    ----------
    qualifier_spec
        A list of tuples.
        Each tuple `t` has 3 string elements.
        `t[0]` specifies the feature type to affect. If this is 'all', all
        feature types are affected.
        `t[1]` gives the key of the qualifier.
        `t[2]` gives the values the qualifier is set to.
    """
    for ftype, k, v in qualifier_spec:
        ftype = ftype.lower()  # not case sensitive
        for feature in record.features:
            # check feature type without case sensitivity
            ftype_have = feature.type.lower() if feature.type else None
            if ftype_have == ftype or ftype == "all":
                # qualifier values should always be lists
                feature.qualifiers[k] = [v]


def add_features(record: SeqRecord, feature_specs: list[tuple[str, str, str, str]]):
    for feature_spec in feature_specs:
        record.features.append(add_feature(record, feature_spec))


def expand_features(record: SeqRecord, expand_specs: list[tuple[str, str]]):
    """
    Modify locations for features in record in place.

    Parameters
    ----------
    record
        SeqRecord whose features are modified.
    expand_specs
        A list of tuples.
        Each tuple `t` has 2 string elements.
        `t[0]` specifies the feature type to affect. If this is 'all', all
        feature types are affected.
        `t[1]` gives the number of bases to expand in either direction.
    """
    # Shrinking features by a large amount may delete them entirely. Therefore,
    # the list of features for the record needs to be rebuilt from scratch,
    # omitting deleted features.
    for feature in record.features:
        for expand_spec in expand_specs:
            ftype, amount = expand_spec
            amount = int(amount)  # arrives from CLI as str
            if feature.type != ftype or feature.location is None:
                continue
            intervals: list[tuple[int, int]] = []
            for part in feature.location.parts:
                start = part.start - amount  # type: ignore
                if start < 0:
                    start = 0
                end = part.end + amount  # type: ignore
                if end > len(record):
                    end = len(record)
                if end >= start:
                    intervals.append((start, end))
            if not intervals:
                # no parts with non-zero size are left after shrinking
                feature.location = None
                # mark for deletion.
                # Do not rely on feature.location in case input features had
                # no location to start with.
                feature.qualifiers["__gbtask_expand_delete"] = True
                continue
            intervals = reduce_intervals(intervals)
            # TODO: factor out intervals2location
            locs = [SimpleLocation(i[0], i[1], feature.strand) for i in intervals]
            if len(locs) > 1:
                feature.location = CompoundLocation(locs)
            else:
                feature.location = locs[0]
    record.features = [
        f for f in record.features if not f.qualifiers.get("__gbtask_expand_delete")
    ]


def _featurepos_str2int(record: SeqRecord, pos: str) -> int:
    """
    Convert positions for feature creation given on the command line as strings
    to integers.
    Negative values are counted from the back of the record.

    Examples
    --------
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> sr = SeqRecord(Seq(100 * "N"))

    >>> _featurepos_str2int(sr, "0")
    0
    >>> _featurepos_str2int(sr, "-0")
    100
    >>> _featurepos_str2int(sr, "1")
    1
    >>> _featurepos_str2int(sr, "-1")
    99
    """
    # set flag to distinguish '0' and '-0'
    if pos[0] == "-":
        negative = True
    else:
        negative = False
    i = abs(int(pos))
    if negative:
        i = len(record) - i
    return i


def add_feature(record: SeqRecord, feature_spec: tuple[str, str, str, str]):
    # convert string CLI arguments to int
    start, end, strand, type = feature_spec
    strand = None if strand == "None" else int(strand)
    start = _featurepos_str2int(record, start)
    end = _featurepos_str2int(record, end)
    if start > end:
        start, end = end, start
    return SeqFeature(SimpleLocation(start, end, strand), type=type)
