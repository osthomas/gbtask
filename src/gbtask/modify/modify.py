from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord


def modify(
    record: SeqRecord,
    pad_left: str,
    pad_right: str,
    include: list[str],
    exclude: list[str],
    features: list[tuple[str, str, str, str]],
    annotations: list[tuple[str, str]],
    qualifiers: list[tuple[str, str, str]],
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
    """
    if not include and not exclude:
        # include everything
        return True
    if include and feature.type in include:
        return True
    if exclude and feature.type not in exclude:
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
        for feature in record.features:
            # check feature type without case sensitivity
            ftype_have = feature.type.lower() if feature.type else None
            if ftype_have == ftype or ftype == "all":
                # qualifier values should always be lists
                feature.qualifiers[k] = [v]


def add_features(record: SeqRecord, feature_specs: list[tuple[str, str, str, str]]):
    for feature_spec in feature_specs:
        record.features.append(add_feature(record, feature_spec))


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
