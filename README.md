# HPA RNA specificity → RDF (hpao)

Human Protein Atlas の `proteinatlas.tsv`（最新リリースのみ対応、現在 v25.1）から、
RNA tissue / single cell type specificity を RDF 化する。

## ファイル

| パス | 内容 |
|---|---|
| `hpa2rdf.py` | 変換スクリプト（Python 3 標準ライブラリのみ） |
| `ontology/hpao.ttl` | オントロジー（クラス・プロパティ・カテゴリ概念） |
| `mapping/tissue.tsv`, `mapping/cell_type.tsv` | HPA 名 → HPA ID（`HPAT_*` / `HPAC_*`）、CALOHA（`caloha_match` 列: exact / close）、UBERON（組織）/ CL（細胞型） |
| `rdf-config/hpa/` | rdf-config（`model.yaml`, `prefix.yaml`, `schema.svg`） |
| `data/proteinatlas.tsv.zip` | 入力（ダウンロードしたもの） |
| `out/hpa.ttl` | 出力 |

## 使い方

```sh
curl -L -o data/proteinatlas.tsv.zip https://www.proteinatlas.org/download/proteinatlas.tsv.zip
python3 hpa2rdf.py data/proteinatlas.tsv.zip -o out/hpa.ttl --hpa-version 25.1
# out/hpa.ttl と ontology/hpao.ttl をトリプルストアに読み込む
```

新しいリリースで分類値や列名が変わった場合はエラーで止まる。
対応表にない組織・細胞型名があるとエラーで止まる。その場合は `mapping/` に新しい番号（既存の最大番号の次）で行を追加する。
対応表には v25.1 に出てこない名前の行も ID 保持のため残してある（出力されるのは使われている用語のみ）。

## データモデル

```turtle
@prefix hpao: <http://med2rdf.org/hpa/ontology#> .
@prefix hpa:  <http://med2rdf.org/hpa/resource/> .

ensembl:ENSG00000000003 a obo:SO_0000704 ;
    hpao:tissueSpecificity hpao:LowTissueSpecificity ;
    hpao:tissueDistribution hpao:DetectedInMany ;
    hpao:cellTypeSpecificity hpao:CellTypeGroupEnriched ;
    hpao:cellTypeDistribution hpao:DetectedInMany ;
    hpao:cellTypeSpecificityScore "5"^^xsd:decimal ;
    hpao:groupEnrichedIn hpao:HPAC_000047 ,
        hpao:HPAC_000070 ;
    hpao:rnaExpression [ hpao:anatomicalEntity hpao:HPAC_000047 ; hpao:nCPM 1280.5 ] ,
        [ hpao:anatomicalEntity hpao:HPAC_000070 ; hpao:nCPM 3513.3 ] .

hpao:HPAC_000070 a hpao:SingleCellType ;
    dcterms:identifier "HPAC_000070" ;
    rdfs:label "Late spermatids" ;
    skos:closeMatch caloha:TS-2535 , obo:CL_4030037 .
```

| 述語 | 値 | 対象 |
|---|---|---|
| `hpao:tissueSpecificity` / `hpao:cellTypeSpecificity` | `TissueEnriched`, `TissueGroupEnriched`, `TissueEnhanced`, `LowTissueSpecificity`, `TissueNotDetected`（細胞型は `CellType…`） | 全遺伝子 |
| `hpao:tissueDistribution` / `hpao:cellTypeDistribution` | `DetectedInAll` / `Many` / `Some` / `Single`, `NotDetected` | 全遺伝子 |
| `hpao:tissueSpecificityScore` / `hpao:cellTypeSpecificityScore` | xsd:decimal（次点の組織・細胞型に対する倍率） | enriched / group enriched のみ |
| `hpao:enrichedIn` / `hpao:groupEnrichedIn` / `hpao:enhancedIn` | `hpao:Tissue` / `hpao:SingleCellType` | specific 列に列挙されたもの |
| `hpao:rnaExpression` → `hpao:anatomicalEntity` + `hpao:nTPM` / `hpao:nCPM` | xsd:decimal | 同上（specific 列の値） |

- 遺伝子の型は `obo:SO_0000704`（gene）
- 組織と細胞型は同じ 3 述語を使い、目的語のクラスで区別する
- 組織・細胞型の IRI は HPA ID（`hpao:HPAT_000018` = liver、`hpao:HPAC_000063` = Hepatocytes）。
  ID と CL / UBERON 対応は参考文献の池田さんのリポジトリの対応表を元にしている（同じ IRI を指す）
- CALOHA は exact / close を区別、UBERON / CL はすべて `skos:closeMatch`
- 組織・細胞型のラベルは HPA の名前そのまま（`skin 1` などの番号付きも含む）
- 収録範囲は `proteinatlas.tsv` にある情報のみ。発現量は specific 列の値なので、
  enriched / group enriched / enhanced の組織・細胞型の分だけ（全組織 × 全遺伝子の行列ではない）

## v25.1 の件数

| | 組織 | 細胞型 |
|---|---|---|
| 遺伝子 | 20,162 | 20,162 |
| enriched / group enriched / enhanced（遺伝子数） | 3,132 / 1,547 / 6,356 | 1,906 / 1,995 / 13,393 |
| Low / Not detected（遺伝子数） | 8,096 / 1,031 | 1,590 / 1,278 |
| 遺伝子–組織リンク（＝発現量の値の数） | 17,870 | 87,130 |
| 用語数 | 37 | 154 |
| CALOHA / UBERON・CL 対応あり | 37 / 37 | 71 / 101（どちらか 107） |

## 参考

- [sh-ikeda/hpa-expr-spec-rdf-converter](https://github.com/sh-ikeda/hpa-expr-spec-rdf-converter) —
  `proteinatlas.tsv` を RDF 化する gawk スクリプト。組織・細胞型の ID と CALOHA / UBERON / CL 対応表
  （`input/hpa_tissue_caloha.tsv`, `input/hpa_cell_cl.tsv`, `input/hpa_tissue_uberon.tsv`）
