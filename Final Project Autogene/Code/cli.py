import argparse
from typing import List

from ecocyc_client import DEFAULT_PGDB_DIR, EcoCycClient, ResolvedGene, get_node_label, resolve_gene_input
from graph_builder import build_graph, export_graph_json


def main() -> None:
    parser = argparse.ArgumentParser(description="Build EcoCyc regulation graph and export JSON for the circuit viewer.")
    parser.add_argument("--pgdb-dir", default=DEFAULT_PGDB_DIR, help="Path to local EcoCyc PGDB directory (default: 29.6).")
    parser.add_argument("-o", "--output-json", default="circuit_data.json", help="Path to write graph JSON consumed by circuit_viewer/index.html (default: circuit_data.json).")
    parser.add_argument("genes", nargs="*", metavar="GENE", help="Gene names or EcoCyc ids; defaults to lacI lacZ.")
    args = parser.parse_args()

    client = EcoCycClient(args.pgdb_dir)
    tokens = list(args.genes) if args.genes else ["LacI", "lacZ"]
    resolved: List[ResolvedGene] = []
    for tok in tokens:
        r = resolve_gene_input(client, tok)
        if r.source == "ecocyc":
            print(f"Resolved {tok!r} -> {r.ecocyc_id}")
        resolved.append(r)

    G, uniprot_labels = build_graph(client, resolved)
    print("\n=== Graph Summary ===")
    print("Nodes:", G.number_of_nodes())
    print("Edges:", G.number_of_edges())

    print("\n=== Edges (Readable) ===")
    name_cache = {}

    def cached_name(nid):
        if nid not in name_cache:
            name_cache[nid] = get_node_label(client, nid, uniprot_labels)
        return name_cache[nid]

    for u, v, d in G.edges(data=True):
        print(f"{cached_name(u)} --{d['interaction']}--> {cached_name(v)}")

    export_graph_json(client, G, args.output_json, name_cache=name_cache, uniprot_labels=uniprot_labels)


if __name__ == "__main__":
    main()
