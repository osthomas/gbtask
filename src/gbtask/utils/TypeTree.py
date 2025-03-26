import networkx as nx

_hierarchy = nx.DiGraph()
_hierarchy.add_node("gene", gtf_level="gene_id")
for rna in ("mRNA", "tRNA", "snRNA", "snoRNA", "lncRNA", "ncRNA"):
    _hierarchy.add_node(rna, gtf_level="transcript_id")
    _hierarchy.add_edge("gene", rna)
    _hierarchy.add_edge(rna, "exon")
_hierarchy.add_node("exon")
_hierarchy.add_node("CDS")
_hierarchy.add_edge("mRNA", "CDS")
_hierarchy.add_node("5'UTR", alias={"gff": "five_prime_utr", "gtf": "UTR"})
_hierarchy.add_edge("mRNA", "5'UTR")
_hierarchy.add_node("3'UTR", alias={"gff": "three_prime_utr", "gtf": "UTR"})
_hierarchy.add_edge("mRNA", "3'UTR")


class TypeTreeError(Exception):
    pass


class TypeTree:
    def __init__(self, hierarchy: nx.DiGraph = _hierarchy):
        self._tree: nx.DiGraph = hierarchy
        self._validate()
        self.root_type: str = self._root_types[0]

    def _validate(self):
        nroots = len(self._root_types)
        if nroots != 1:
            msg = f"The type hierarchy can only have one root, but it has{nroots} ({self._root_types})"
            raise TypeTreeError(msg)

    @property
    def _root_types(self) -> list[str]:
        """
        Return the types of root nodes of the underlying feature type tree.
        """
        return list([node for node, deg in self._tree.in_degree() if deg == 0])

    @property
    def root_nodes(self) -> list:
        """
        Return root nodes of the underlying feature type tree.
        """
        return [self._tree[rt] for rt in self.root_types]

    def parent_types(self, feature_type: str):
        try:
            return self._tree.predecessors(feature_type)
        except nx.NetworkXError:
            # requested feature type is not in the type hierarchy
            return tuple()

    def alias(self, feature_type: str, format: str):
        """
        The alias definitions are fetched from the nodes' metadata (key 'alias').
        This should be a dict[str, str], with keys corresponding to file format names
        (eg. "gff", "gtf").
        If no aliases are available for a node, or if it there is no entry for the
        requested format type, the original label is retained.
        """
        try:
            alias = self._tree.nodes[feature_type]["alias"][format]
        except KeyError:
            alias = feature_type  # retain original label by default
        return alias

    def relabel(self, format: str):
        """
        Relabel node labels to conform to feature type names of different file formats.

        The alias definitions are fetched from the nodes' metadata (key 'alias').
        This should be a dict[str, str], with keys corresponding to file format names
        (eg. "gff", "gtf").
        If no aliases are available for a node, or if it there is no entry for the
        requested format type, the original label is retained.
        """
        mapping = {}
        for orig in self._tree.nodes:
            alias = self.alias(orig, format)
            if alias in mapping.values():
                raise NotImplementedError(
                    f"Type name alias {alias} for {orig} is not unique!"
                )
            else:
                mapping[orig] = alias
        return self.__class__(nx.relabel_nodes(self._tree, mapping))

    @property
    def gtf_levels(self):
        """
        Return GTF grouping levels available in the tree.
        """
        levels = []
        for n in self._tree:
            try:
                level = self._tree.nodes[n]["gtf_level"]
                if level not in levels:
                    levels.append(level)
            except KeyError:
                pass
        return levels

    def gtf_level_types(self, level: str):
        """
        Return feature types at which `level` should be set for the respective
        feature itself and its descendants.

        Parameters

        level
            GTF level to fetch types for, typically "gene_id", "transcript_id"
        """
        types = []
        for node in self._tree:
            try:
                if self._tree.nodes[node]["gtf_level"] == level and node not in types:
                    types.append(node)
            except KeyError:
                pass
        return types


DefaultTypeTree = TypeTree()


def plot(tree: TypeTree, ax=None):
    """
    Plot a TypeTree's underlying feature tree.
    """
    from matplotlib import pyplot as plt

    g = tree._tree
    pos = nx.nx_agraph.graphviz_layout(g, prog="dot")
    nx.draw_networkx_edges(g, pos, ax=ax, edge_color="black", node_size=500)
    nx.draw_networkx_labels(g, pos, ax=ax, font_size=10, bbox={"ec": "k", "fc": "w"})
    return plt.gca()
