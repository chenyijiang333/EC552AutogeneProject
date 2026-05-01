"""Backward-compatible facade for the refactored Autogene modules."""

from cli import main
from ecocyc_client import (
    DEFAULT_PGDB_DIR,
    ORG_ID,
    EcoCycClient,
    ResolvedGene,
    ecocyc_node_uniprot_entry_url,
    get_node_label,
    resolve_gene_input,
)
from graph_builder import (
    _node_kind,
    build_graph,
    export_graph_json,
    seed_graph_node_id,
)
from remote_clients import (
    GENE_NODE_PREFIX,
    UNIPROT_NODE_PREFIX,
    normalize_legacy_uniprot_gene_search_url,
    uniprot_entry_url,
    uniprot_url_for_graph_node,
)
from sequence_analysis import analyze_sequences_for_graph_hints


if __name__ == "__main__":
    main()
