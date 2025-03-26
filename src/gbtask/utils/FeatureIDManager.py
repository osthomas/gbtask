from collections import Counter

from Bio.SeqFeature import SeqFeature


class FeatureIDManager:
    def __init__(self, features: list[SeqFeature], prefix: str = ""):
        """
        Manage IDs for groups of features, ensuring disambiguation of non-unique
        feature IDs.


        Parameters
        ----------

        features
            A list of SeqFeatures for which to generate IDs.
        prefix
            Prefix to prepend to all generated IDs. Intended for marking group
            membership of features associated with the FeatureIDManager.
        """
        # Sort by start position for successive numbering
        features = sorted(features, key=lambda f: f.location.start if f.location else 0)
        # Do not label the *same* feature repeatedly
        # Same means stored at the same physical address, not equal as Python
        # object.
        # Distinct feature objects compare as equal if their attributes are
        # compatible. This is a problem for newly instantiated features that
        # are indistinguishable from other features.
        seen = []
        self.features: list[SeqFeature] = []
        for feature in features:
            if id(feature) not in seen:
                seen.append(id(feature))
                self.features.append(feature)
        self.prefix: str = prefix
        self._counter: Counter = Counter()

    def set_feature_ids(self):
        self._counter = Counter()  # reset counter
        for feature in self.features:
            self._set_feature_id(feature)

    def _set_feature_id(
        self,
        feature: SeqFeature,
    ):
        """
        Set feature IDs to `feature.id` and `feature.qualifiers["ID"]`.

        Parameters
        ----------

        feature
            This SeqFeature is modified in place!

        Return
        ------
        The new feature ID
        """
        base = feature.type or "Unknown"
        self._counter[base] += 1
        if self.prefix:
            final = f"{self.prefix}:{base}"
        else:
            final = base
        # Disambiguate if not unique for whatever reason
        final = f"{final}.{self._counter[base]:02}"
        feature.id = final
        feature.qualifiers["ID"] = [final]
        return final
