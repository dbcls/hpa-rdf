#!/usr/bin/env python3
"""Convert Human Protein Atlas proteinatlas.tsv (latest release) into RDF.

Vocabulary: hpao.ttl (http://med2rdf.org/hpao/ontology#)

  ensembl:ENSG... a obo:SO_0000704 ;
      hpao:tissueSpecificity hpao:TissueEnriched ;
      hpao:tissueDistribution hpao:DetectedInSingle ;
      hpao:tissueSpecificityScore 16 ;
      hpao:enrichedIn hpat:tissue_liver ;
      hpao:rnaExpression [ hpao:anatomicalEntity hpat:tissue_liver ; hpao:nTPM 1234.5 ] ;
      hpao:cellTypeSpecificity hpao:CellTypeEnhanced ;
      ...
  hpat:tissue_liver a hpao:Tissue ; rdfs:label "liver" ; skos:exactMatch caloha:TS-0564 .

Every gene row gets a specificity and distribution category (incl. Low / Not
detected). enrichedIn / groupEnrichedIn / enhancedIn link to the tissues / cell
types listed in "RNA tissue specific nTPM" / "RNA single cell type specific nCPM",
and hpao:rnaExpression holds the value given there.

Usage:
  hpa2rdf.py proteinatlas.tsv[.zip] -o hpa.ttl [--hpa-version 25.1]
"""

import argparse
import csv
import datetime
import io
import re
import sys
import unicodedata
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent

PREFIXES = """@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .
@prefix skos: <http://www.w3.org/2004/02/skos/core#> .
@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix pav: <http://purl.org/pav/> .
@prefix hpao: <http://med2rdf.org/hpao/ontology#> .
@prefix hpat: <http://med2rdf.org/hpao/term/> .
@prefix hpa: <http://med2rdf.org/hpao/resource/> .
@prefix ensembl: <http://identifiers.org/ensembl/> .
@prefix caloha: <http://purl.obolibrary.org/obo/caloha.obo#> .
@prefix obo: <http://purl.obolibrary.org/obo/> .
"""

GROUP_PREDICATE = {"enriched": "hpao:enrichedIn", "group": "hpao:groupEnrichedIn", "enhanced": "hpao:enhancedIn"}

DISTRIBUTION = {
    "Detected in all": "hpao:DetectedInAll",
    "Detected in many": "hpao:DetectedInMany",
    "Detected in some": "hpao:DetectedInSome",
    "Detected in single": "hpao:DetectedInSingle",
    "Not detected": "hpao:NotDetected",
}

# One entry per HPA classification.
#   category value in proteinatlas.tsv -> (category concept, predicate group or None)
DATASETS = [
    {
        "column": "RNA tissue",
        "values": "RNA tissue specific nTPM",
        "class": "hpao:Tissue",
        "term_prefix": "tissue_",
        "mapping": "tissue_caloha.tsv",
        "specificity": "hpao:tissueSpecificity",
        "distribution": "hpao:tissueDistribution",
        "score": "hpao:tissueSpecificityScore",
        "unit": "hpao:nTPM",
        "categories": {
            "Tissue enriched": ("hpao:TissueEnriched", "enriched"),
            "Group enriched": ("hpao:TissueGroupEnriched", "group"),
            "Tissue enhanced": ("hpao:TissueEnhanced", "enhanced"),
            "Low tissue specificity": ("hpao:LowTissueSpecificity", None),
            "Not detected": ("hpao:TissueNotDetected", None),
        },
    },
    {
        "column": "RNA single cell type",
        "values": "RNA single cell type specific nCPM",
        "class": "hpao:SingleCellType",
        "term_prefix": "celltype_",
        "mapping": "cell_type_caloha.tsv",
        "specificity": "hpao:cellTypeSpecificity",
        "distribution": "hpao:cellTypeDistribution",
        "score": "hpao:cellTypeSpecificityScore",
        "unit": "hpao:nCPM",
        "categories": {
            "Cell type enriched": ("hpao:CellTypeEnriched", "enriched"),
            "Group enriched": ("hpao:CellTypeGroupEnriched", "group"),
            "Cell type enhanced": ("hpao:CellTypeEnhanced", "enhanced"),
            "Low cell type specificity": ("hpao:LowCellTypeSpecificity", None),
            "Not detected": ("hpao:CellTypeNotDetected", None),
        },
    },
]


def open_tsv(path):
    path = Path(path)
    if path.suffix == ".zip":
        zf = zipfile.ZipFile(path)
        name = next(n for n in zf.namelist() if n.endswith(".tsv"))
        return io.TextIOWrapper(zf.open(name), encoding="utf-8")
    return open(path, encoding="utf-8")


def load_mapping(path):
    """hpa_name -> list of (property, CURIE)"""
    mapping = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            name, caloha_id, _label, match = (line.rstrip("\n").split("\t") + [""] * 4)[:4]
            links = []
            if caloha_id:
                prop = "skos:closeMatch" if match == "close" else "skos:exactMatch"
                links.append((prop, f"caloha:{caloha_id}"))
            mapping[name] = links
    return mapping


def parse_values(value, gene):
    # "Early spermatids: 313.3;Late spermatids: 1149.4" -> [("Early spermatids", "313.3"), ("Late spermatids", "1149.4")]
    if not value:
        return []
    pairs = []
    for kv in value.split(";"):
        name, num = (x.strip() for x in kv.rsplit(":", 1))
        try:
            float(num)
        except ValueError:
            sys.exit(f"{gene}: non-numeric value in {kv!r}")
        pairs.append((name, num if "." in num else num + ".0"))  # keep xsd:decimal
    return pairs


def local_name(name):
    # "Müller glia" -> "Muller_glia", "Pituicytes/FSCs" -> "Pituicytes_FSCs"
    s = unicodedata.normalize("NFKD", name).encode("ascii", "ignore").decode()
    return re.sub(r"[^A-Za-z0-9-]+", "_", s).strip("_")


def literal(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def check(value, allowed, what, gene):
    if value not in allowed:
        sys.exit(f"{gene}: unexpected {what} {value!r}")
    return allowed[value]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="proteinatlas.tsv or proteinatlas.tsv.zip")
    ap.add_argument("-o", "--output", default="hpa.ttl", help="output Turtle file (default: hpa.ttl)")
    ap.add_argument("-m", "--mapping-dir", default=HERE / "mapping", type=Path,
                    help="directory holding tissue_caloha.tsv and cell_type_caloha.tsv")
    ap.add_argument("--hpa-version", help="HPA release, e.g. 25.1 (recorded in the dataset description)")
    args = ap.parse_args()

    with open_tsv(args.input) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        missing = [c for d in DATASETS for c in (d["values"], d["column"] + " specificity") if c not in reader.fieldnames]
        if missing:
            sys.exit(f"columns not found (not the expected HPA release?): {missing}")
    print(f"{args.input}: {len(rows)} genes", file=sys.stderr)

    genes = {}   # ENSG -> list of "predicate object" strings
    terms = {}   # CURIE -> (class, label, links)
    for d in DATASETS:
        mapping = load_mapping(args.mapping_dir / d["mapping"])
        unlisted = set()
        n_links = 0
        for row in rows:
            gene = row["Ensembl"]
            po = genes.setdefault(gene, ["a obo:SO_0000704"])
            category, group = check(row[d["column"] + " specificity"], d["categories"], "specificity", gene)
            po.append(f"{d['specificity']} {category}")
            po.append(f"{d['distribution']} " + check(row[d["column"] + " distribution"], DISTRIBUTION, "distribution", gene))
            score = row[d["column"] + " specificity score"]
            if score:
                po.append(f'{d["score"]} "{float(score):g}"^^xsd:decimal')
            pairs = parse_values(row[d["values"]], gene)
            if pairs and not group:
                sys.exit(f"{gene}: {category} but {d['values']} = {row[d['values']]!r}")
            objs = []
            exprs = []
            for name, num in pairs:
                curie = "hpat:" + d["term_prefix"] + local_name(name)
                prev = terms.setdefault(curie, (d["class"], name, mapping.get(name, [])))
                if prev[1] != name:
                    sys.exit(f"term ID collision: {name!r} and {prev[1]!r} -> {curie}")
                if name not in mapping:
                    unlisted.add(name)
                objs.append(curie)
                exprs.append(f"[ hpao:anatomicalEntity {curie} ; {d['unit']} {num} ]")
            if objs:
                po.append(f"{GROUP_PREDICATE[group]} " + " ,\n        ".join(objs))
                po.append("hpao:rnaExpression " + " ,\n        ".join(exprs))
                n_links += len(objs)
        n_terms = sum(1 for t in terms.values() if t[0] == d["class"])
        n_mapped = sum(1 for t in terms.values() if t[0] == d["class"] and t[2])
        print(f"  {d['column']}: {n_links} links to {n_terms} terms ({n_mapped} with CALOHA match)", file=sys.stderr)
        for name in sorted(unlisted):
            print(f"    not in {d['mapping']}: {name}", file=sys.stderr)

    with open(args.output, "w", encoding="utf-8") as out:
        out.write(PREFIXES + "\n")
        out.write("hpa:dataset a dcterms:Dataset ;\n"
                  '    dcterms:title "Human Protein Atlas RNA tissue and single cell type specificity" ;\n'
                  "    dcterms:source <https://www.proteinatlas.org/about/download> ;\n"
                  "    dcterms:license <https://creativecommons.org/licenses/by/4.0/> ;\n")
        if args.hpa_version:
            out.write(f"    pav:version {literal(args.hpa_version)} ;\n")
        out.write(f'    dcterms:created "{datetime.date.today().isoformat()}"^^xsd:date .\n\n')

        for curie, (cls, label, links) in sorted(terms.items()):
            out.write(f"{curie} a {cls} ;\n    rdfs:label {literal(label)}")
            for prop, target in links:
                out.write(f" ;\n    {prop} {target}")
            out.write(" .\n")
        out.write("\n")

        for gene, po in genes.items():
            out.write(f"ensembl:{gene} " + " ;\n    ".join(po) + " .\n")

    print(f"wrote {args.output}: {len(genes)} genes, {len(terms)} terms", file=sys.stderr)


if __name__ == "__main__":
    main()
