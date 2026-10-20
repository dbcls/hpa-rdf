#!/usr/bin/env python3
"""Convert Human Protein Atlas RNA specificity (latest release) into RDF.

Vocabulary: ontology/hpao.ttl (http://med2rdf.org/hpa/ontology#)

  ensembl:ENSG... a obo:SO_0000704 ;
      hpao:tissueSpecificity hpao:TissueEnriched ;
      hpao:tissueDistribution hpao:DetectedInSingle ;
      hpao:tissueSpecificityScore 16 ;
      hpao:enrichedIn hpao:HPAT_000018 ;
      hpao:rnaExpression [ hpao:anatomicalEntity hpao:HPAT_000018 ; hpao:nTPM 1234.5 ] ;
      hpao:cellTypeSpecificity hpao:CellTypeEnhanced ;
      hpao:cellTypeGroupSpecificity hpao:CellTypeEnriched ;
      ...
  hpao:HPAC_000063 a owl:Class , hpao:SingleCellType ; dcterms:identifier "HPAC_000063" ;
      rdfs:label "Hepatocytes" ; rdfs:subClassOf hpao:HPACG_000022 ;
      skos:exactMatch caloha:TS-0454 ; skos:closeMatch obo:CL_0000182 .
  hpao:HPACG_000022 a owl:Class , hpao:SingleCellTypeGroup ; rdfs:subClassOf hpao:HPACC_000012 .
  hpao:HPACC_000012 a owl:Class , hpao:SingleCellTypeClass .

Every gene row gets a specificity and distribution category (incl. Low / Not
detected) for tissue, single cell type and single cell type group.
enrichedIn / coEnrichedIn / enhancedIn link to the tissues / cell types / cell
type groups listed in the "... specific nTPM|nCPM" columns, and
hpao:rnaExpression holds the value given there. HPA's "Group enriched" is
hpao:coEnrichedIn, so that "group" only means a cell type group.

Inputs:
  proteinatlas.tsv                      classification and values
  rna_single_cell_type_cell_types.tsv   cell type -> cell type group -> cell type class
  mapping/*.tsv                         HPA IDs (HPAT_ / HPAC_ / HPACG_ / HPACC_) and
                                        CALOHA / UBERON / CL matches

Usage:
  hpa2rdf.py proteinatlas.tsv[.zip] rna_single_cell_type_cell_types.tsv[.zip] -o hpa.ttl [--hpa-version 25.1 --hpa-date 2026-05-25]
"""

import argparse
import csv
import datetime
import io
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent

PREFIXES = """@prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#> .
@prefix owl: <http://www.w3.org/2002/07/owl#> .
@prefix xsd: <http://www.w3.org/2001/XMLSchema#> .
@prefix skos: <http://www.w3.org/2004/02/skos/core#> .
@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix pav: <http://purl.org/pav/> .
@prefix hpao: <http://med2rdf.org/hpa/ontology#> .
@prefix hpa: <http://med2rdf.org/hpa/resource/> .
@prefix ensembl: <http://identifiers.org/ensembl/> .
@prefix caloha: <http://purl.obolibrary.org/obo/caloha.obo#> .
@prefix obo: <http://purl.obolibrary.org/obo/> .
"""

LINK_PREDICATE = {"enriched": "hpao:enrichedIn", "co": "hpao:coEnrichedIn", "enhanced": "hpao:enhancedIn"}

DISTRIBUTION = {
    "Detected in all": "hpao:DetectedInAll",
    "Detected in many": "hpao:DetectedInMany",
    "Detected in some": "hpao:DetectedInSome",
    "Detected in single": "hpao:DetectedInSingle",
    "Not detected": "hpao:NotDetected",
}

CELL_TYPE_CATEGORIES = {
    "Cell type enriched": ("hpao:CellTypeEnriched", "enriched"),
    "Group enriched": ("hpao:CellTypeCoEnriched", "co"),
    "Cell type enhanced": ("hpao:CellTypeEnhanced", "enhanced"),
    "Low cell type specificity": ("hpao:LowCellTypeSpecificity", None),
    "Not detected": ("hpao:CellTypeNotDetected", None),
}

# Anatomical levels: mapping file -> class of its terms
LEVELS = {
    "tissue": ("tissue.tsv", "hpao:Tissue"),
    "cell_type": ("cell_type.tsv", "hpao:SingleCellType"),
    "group": ("cell_type_group.tsv", "hpao:SingleCellTypeGroup"),
    "class": ("cell_type_class.tsv", "hpao:SingleCellTypeClass"),
}

# One entry per HPA classification in proteinatlas.tsv.
#   category value -> (category concept, link predicate key or None)
DATASETS = [
    {
        "column": "RNA tissue",
        "values": "RNA tissue specific nTPM",
        "level": "tissue",
        "specificity": "hpao:tissueSpecificity",
        "distribution": "hpao:tissueDistribution",
        "score": "hpao:tissueSpecificityScore",
        "unit": "hpao:nTPM",
        "categories": {
            "Tissue enriched": ("hpao:TissueEnriched", "enriched"),
            "Group enriched": ("hpao:TissueCoEnriched", "co"),
            "Tissue enhanced": ("hpao:TissueEnhanced", "enhanced"),
            "Low tissue specificity": ("hpao:LowTissueSpecificity", None),
            "Not detected": ("hpao:TissueNotDetected", None),
        },
    },
    {
        "column": "RNA single cell type",
        "values": "RNA single cell type specific nCPM",
        "level": "cell_type",
        "specificity": "hpao:cellTypeSpecificity",
        "distribution": "hpao:cellTypeDistribution",
        "score": "hpao:cellTypeSpecificityScore",
        "unit": "hpao:nCPM",
        "categories": CELL_TYPE_CATEGORIES,
    },
    {
        "column": "RNA single cell type group",
        "values": "RNA single cell type group specific nCPM",
        "level": "group",
        "specificity": "hpao:cellTypeGroupSpecificity",
        "distribution": "hpao:cellTypeGroupDistribution",
        "score": "hpao:cellTypeGroupSpecificityScore",
        "unit": "hpao:nCPM",
        "categories": CELL_TYPE_CATEGORIES,  # HPA uses the cell type categories at group level
    },
]


def open_tsv(path):
    path = Path(path)
    if path.suffix == ".zip":
        zf = zipfile.ZipFile(path)
        name = next(n for n in zf.namelist() if n.endswith(".tsv"))
        return io.TextIOWrapper(zf.open(name), encoding="utf-8")
    return open(path, encoding="utf-8")


def load_mapping(path, cls):
    """lower-cased hpa_name -> term dict

    Columns are read by header name: hpa_id, hpa_name, and optionally
    caloha_id + caloha_match (exact / close) and other *_id columns (UBERON / CL, closeMatch).
    """
    terms = {}
    with open(path, encoding="utf-8") as f:
        header = f.readline().lstrip("# ").rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t") + [""] * len(header)))
            links = []
            if r.get("caloha_id"):
                prop = "skos:closeMatch" if r.get("caloha_match") == "close" else "skos:exactMatch"
                links.append((prop, f"caloha:{r['caloha_id']}"))
            for col in header:
                if col.endswith("_id") and col not in ("hpa_id", "caloha_id") and r[col]:
                    links.append(("skos:closeMatch", f"obo:{r[col]}"))
            terms[r["hpa_name"].lower()] = {"curie": "hpao:" + r["hpa_id"], "id": r["hpa_id"], "label": r["hpa_name"],
                                            "class": cls, "links": links, "parent": None}
    return terms


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


def release_date(s):
    try:
        return datetime.datetime.strptime(s, "%Y-%m-%d").date()
    except ValueError:
        raise argparse.ArgumentTypeError(f"not YYYY-MM-DD: {s!r}")


def literal(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def check(value, allowed, what, gene):
    if value not in allowed:
        sys.exit(f"{gene}: unexpected {what} {value!r}")
    return allowed[value]


def lookup(terms, name, level, where):
    term = terms[level].get(name.lower())
    if not term:
        sys.exit(f"{name!r} ({where}) not in mapping/{LEVELS[level][0]} (add a row with a new HPA ID)")
    return term


def load_hierarchy(path, terms):
    """Set cell type -> group -> class parents from rna_single_cell_type_cell_types.tsv."""
    with open_tsv(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            cell = lookup(terms, r["Cell type"], "cell_type", path)
            group = lookup(terms, r["Cell type group"], "group", path)
            klass = lookup(terms, r["Cell type class"], "class", path)
            for child, parent in ((cell, group), (group, klass)):
                if child["parent"] not in (None, parent["curie"]):
                    sys.exit(f"{child['label']!r} has two parents: {child['parent']} and {parent['curie']}")
                child["parent"] = parent["curie"]
    missing = [t["label"] for level in ("cell_type", "group") for t in terms[level].values() if not t["parent"]]
    if missing:
        sys.exit(f"not in {path}: {', '.join(missing)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="proteinatlas.tsv or proteinatlas.tsv.zip")
    ap.add_argument("cell_types", help="rna_single_cell_type_cell_types.tsv[.zip] (cell type hierarchy)")
    ap.add_argument("-o", "--output", default="hpa.ttl", help="output Turtle file (default: hpa.ttl)")
    ap.add_argument("-m", "--mapping-dir", default=HERE / "mapping", type=Path,
                    help="directory holding " + ", ".join(f for f, _ in LEVELS.values()))
    ap.add_argument("--hpa-version", help="HPA release, e.g. 25.1 (recorded in the dataset description)")
    ap.add_argument("--hpa-date", type=release_date,
                    help="HPA release date, e.g. 2026-05-25 (recorded as dcterms:issued)")
    args = ap.parse_args()

    terms = {level: load_mapping(args.mapping_dir / f, cls) for level, (f, cls) in LEVELS.items()}
    load_hierarchy(args.cell_types, terms)

    with open_tsv(args.input) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        missing = [c for d in DATASETS for c in (d["values"], d["column"] + " specificity") if c not in reader.fieldnames]
        if missing:
            sys.exit(f"columns not found (not the expected HPA release?): {missing}")
    print(f"{args.input}: {len(rows)} genes", file=sys.stderr)

    genes = {}   # ENSG -> list of "predicate object" strings
    used = set()
    for d in DATASETS:
        n_links = 0
        for row in rows:
            gene = row["Ensembl"]
            po = genes.setdefault(gene, ["a obo:SO_0000704"])
            category, link = check(row[d["column"] + " specificity"], d["categories"], "specificity", gene)
            po.append(f"{d['specificity']} {category}")
            po.append(f"{d['distribution']} " + check(row[d["column"] + " distribution"], DISTRIBUTION, "distribution", gene))
            score = row[d["column"] + " specificity score"]
            if score:
                po.append(f'{d["score"]} "{float(score):g}"^^xsd:decimal')
            pairs = parse_values(row[d["values"]], gene)
            if pairs and not link:
                sys.exit(f"{gene}: {category} but {d['values']} = {row[d['values']]!r}")
            objs = []
            exprs = []
            for name, num in pairs:
                curie = lookup(terms, name, d["level"], d["values"])["curie"]
                used.add(curie)
                objs.append(curie)
                exprs.append(f"[ hpao:anatomicalEntity {curie} ; {d['unit']} {num} ]")
            if objs:
                po.append(f"{LINK_PREDICATE[link]} " + " ,\n        ".join(objs))
                po.append("hpao:rnaExpression " + " ,\n        ".join(exprs))
                n_links += len(objs)
        print(f"  {d['column']}: {n_links} links", file=sys.stderr)

    # tissues: those used; cell types / groups / classes: the whole hierarchy
    out_terms = [t for t in terms["tissue"].values() if t["curie"] in used]
    out_terms += [t for level in ("cell_type", "group", "class") for t in terms[level].values()
                  if t["curie"] in used or t["parent"] or level == "class"]
    for level, (_, cls) in LEVELS.items():
        ts = [t for t in out_terms if t["class"] == cls]
        print(f"  {cls}: {len(ts)} terms ({sum(1 for t in ts if t['links'])} with external match)", file=sys.stderr)

    with open(args.output, "w", encoding="utf-8") as out:
        out.write(PREFIXES + "\n")
        out.write("hpa:dataset a dcterms:Dataset ;\n"
                  '    dcterms:title "Human Protein Atlas RNA tissue, single cell type and cell type group specificity" ;\n'
                  "    dcterms:source <https://www.proteinatlas.org/about/download> ;\n"
                  "    dcterms:license <https://creativecommons.org/licenses/by/4.0/> ;\n")
        if args.hpa_version:
            out.write(f"    pav:version {literal(args.hpa_version)} ;\n")
        if args.hpa_date:
            out.write(f'    dcterms:issued "{args.hpa_date.isoformat()}"^^xsd:date ;\n')
        out.write(f'    dcterms:created "{datetime.date.today().isoformat()}"^^xsd:date .\n\n')

        for t in sorted(out_terms, key=lambda t: t["id"]):
            out.write(f"{t['curie']} a owl:Class , {t['class']} ;\n"
                      f"    dcterms:identifier {literal(t['id'])} ;\n    rdfs:label {literal(t['label'])}")
            if t["parent"]:
                out.write(f" ;\n    rdfs:subClassOf {t['parent']}")
            for prop, target in t["links"]:
                out.write(f" ;\n    {prop} {target}")
            out.write(" .\n")
        out.write("\n")

        for gene, po in genes.items():
            out.write(f"ensembl:{gene} " + " ;\n    ".join(po) + " .\n")

    print(f"wrote {args.output}: {len(genes)} genes, {len(out_terms)} terms", file=sys.stderr)


if __name__ == "__main__":
    main()
