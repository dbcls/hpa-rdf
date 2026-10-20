# HPA RNA specificity → RDF (hpao)

Human Protein Atlas（最新リリースのみ対応、現在 v25.1）の RNA tissue / single cell type /
single cell type group specificity を RDF 化する。

## ファイル

| パス | 内容 |
|---|---|
| `hpa2rdf.py` | 変換スクリプト（Python 3 標準ライブラリのみ） |
| `ontology/hpao.ttl` | オントロジー（クラス・プロパティ・カテゴリ概念） |
| `mapping/tissue.tsv` | 組織 → `HPAT_*`、CALOHA（`caloha_match` 列: exact / close）、UBERON |
| `mapping/cell_type.tsv` | 細胞型 → `HPAC_*`、CALOHA、CL |
| `mapping/cell_type_group.tsv` | 細胞型グループ → `HPACG_*`、CL |
| `mapping/cell_type_class.tsv` | 細胞型クラス → `HPACC_*` |
| `rdf-config/hpa/` | rdf-config（`model.yaml`, `prefix.yaml`, `schema.svg`） |
| `data/` | 入力（ダウンロードしたもの） |
| `out/hpa.ttl` | 出力 |

## 使い方

```sh
curl -L -o data/proteinatlas.tsv.zip https://www.proteinatlas.org/download/proteinatlas.tsv.zip
curl -L -o data/rna_single_cell_type_cell_types.tsv.zip \
    https://www.proteinatlas.org/download/tsv/rna_single_cell_type_cell_types.tsv.zip
python3 hpa2rdf.py data/proteinatlas.tsv.zip data/rna_single_cell_type_cell_types.tsv.zip \
    -o out/hpa.ttl --hpa-version 25.1 --hpa-date 2026-05-25
# out/hpa.ttl と ontology/hpao.ttl をトリプルストアに読み込む
```

| 入力 | 使う内容 |
|---|---|
| `proteinatlas.tsv` | 組織・細胞型・細胞型グループの specificity / distribution / specificity score / specific nTPM・nCPM |
| `rna_single_cell_type_cell_types.tsv` | 細胞型 → 細胞型グループ → 細胞型クラス（[Single cell type data](https://www.proteinatlas.org/humanproteome/single+cell/single+cell+type/data)） |

`--hpa-version` と `--hpa-date` は HPA のリリース番号と公開日で、データセット記述の `pav:version` と `dcterms:issued` になる
（`dcterms:created` は変換した日）。
新しいリリースで分類値や列名が変わった場合はエラーで止まる。
対応表にない名前があるとエラーで止まるので、`mapping/` の該当ファイルに新しい番号（既存の最大番号の次）で行を追加する。
名前の照合は大文字小文字を区別しない（`rna_single_cell_type_cell_types.tsv` は小文字表記）。

## データモデル

```turtle
@prefix hpao: <http://med2rdf.org/hpa/ontology#> .
@prefix hpa:  <http://med2rdf.org/hpa/resource/> .

ensembl:ENSG00000000003 a obo:SO_0000704 ;
    hpao:tissueSpecificity hpao:LowSpecificity ;
    hpao:tissueDistribution hpao:DetectedInMany ;
    hpao:cellTypeSpecificity hpao:CoEnriched ;
    hpao:cellTypeDistribution hpao:DetectedInMany ;
    hpao:cellTypeSpecificityScore "5"^^xsd:decimal ;
    hpao:coEnrichedIn hpao:HPAC_000047 , hpao:HPAC_000070 ;
    hpao:rnaExpression [ hpao:anatomicalEntity hpao:HPAC_000047 ; hpao:nCPM 1280.5 ] ,
        [ hpao:anatomicalEntity hpao:HPAC_000070 ; hpao:nCPM 3513.3 ] ;
    hpao:cellTypeGroupSpecificity hpao:CoEnriched ;
    hpao:cellTypeGroupDistribution hpao:DetectedInMany ;
    hpao:cellTypeGroupSpecificityScore "8"^^xsd:decimal ;
    hpao:coEnrichedIn hpao:HPACG_000005 , hpao:HPACG_000047 ;
    hpao:rnaExpression [ hpao:anatomicalEntity hpao:HPACG_000005 ; hpao:nCPM 1280.5 ] ,
        [ hpao:anatomicalEntity hpao:HPACG_000047 ; hpao:nCPM 3513.3 ] .

hpao:HPAC_000063 a owl:Class , hpao:SingleCellType ;        # 細胞型
    dcterms:identifier "HPAC_000063" ;
    rdfs:label "Hepatocytes" ;
    rdfs:subClassOf hpao:HPACG_000022 ;
    skos:exactMatch caloha:TS-0454 ;
    skos:closeMatch obo:CL_0000182 .

hpao:HPACG_000022 a owl:Class , hpao:SingleCellTypeGroup ;  # 細胞型グループ
    dcterms:identifier "HPACG_000022" ;
    rdfs:label "Hepatocytes" ;
    rdfs:subClassOf hpao:HPACC_000012 ;
    skos:closeMatch obo:CL_0000182 .

hpao:HPACC_000012 a owl:Class , hpao:SingleCellTypeClass ;  # 細胞型クラス
    dcterms:identifier "HPACC_000012" ;
    rdfs:label "specialized epithelial cells" .
```

| 述語 | 値 | 対象 |
|---|---|---|
| `hpao:tissueSpecificity` / `hpao:cellTypeSpecificity` / `hpao:cellTypeGroupSpecificity` | `Enriched`, `CoEnriched`, `Enhanced`, `LowSpecificity`, `NotDetected` | 全遺伝子 |
| `hpao:tissueDistribution` / `hpao:cellTypeDistribution` / `hpao:cellTypeGroupDistribution` | `DetectedInAll` / `Many` / `Some` / `Single`, `NotDetected` | 全遺伝子 |
| `hpao:tissueSpecificityScore` / `hpao:cellTypeSpecificityScore` / `hpao:cellTypeGroupSpecificityScore` | xsd:decimal（次点に対する倍率） | enriched / co-enriched のみ |
| `hpao:enrichedIn` / `hpao:coEnrichedIn` / `hpao:enhancedIn` | `hpao:Tissue` / `hpao:SingleCellType` / `hpao:SingleCellTypeGroup` | specific 列に列挙されたもの |
| `hpao:rnaExpression` → `hpao:anatomicalEntity` + `hpao:nTPM` / `hpao:nCPM` | xsd:decimal | 同上（specific 列の値） |

- 遺伝子の型は `obo:SO_0000704`（gene）
- 組織・細胞型・細胞型グループは同じ 3 述語を使い、目的語の型でどの分類かを区別する
- HPA の "Group enriched"（2〜5 個で 4 倍以上）は `coEnriched`。"group" は細胞型グループの意味にだけ使う
  （`skos:prefLabel` は HPA の表記 "Group enriched" のまま）
- 分類値・分布値は 3 段で共通。どの段の分類かは述語で分かる。HPA の表記（"Tissue enriched" など）は `skos:altLabel`
- `hpao:NotDetected` は分類値と分布値を兼ねる（v25.1 では分類が Not detected ⇔ 分布が Not detected）
- 組織・細胞型・グループ・クラスは `owl:Class` で、段を表す型（`hpao:Tissue` など）も持つ（パニング）。
  細胞型 `rdfs:subClassOf` グループ `rdfs:subClassOf` クラス
- IRI は HPA ID: 組織 `HPAT_`、細胞型 `HPAC_`、グループ `HPACG_`、クラス `HPACC_`（グループ・クラスは名前のアルファベット順）。
  細胞型とグループ、グループとクラスで同じ名前のものがあるので段ごとに別 ID にしている
- 外部対応: CALOHA は exact / close を区別、UBERON / CL はすべて `skos:closeMatch`。クラスには対応なし
- ラベルは HPA の名前そのまま（`skin 1` などの番号付きも含む。クラスは小文字表記）
- 発現量は specific 列の値なので、enriched / co-enriched / enhanced の組織・細胞型・グループの分だけ
  （全組織 × 全遺伝子の行列ではない）

## v25.1 の件数

| | 組織 | 細胞型 | 細胞型グループ |
|---|---|---|---|
| 遺伝子 | 20,162 | 20,162 | 20,162 |
| enriched / co-enriched / enhanced（遺伝子数） | 3,132 / 1,547 / 6,356 | 1,906 / 1,995 / 13,393 | 3,321 / 2,174 / 10,011 |
| Low / Not detected（遺伝子数） | 8,096 / 1,031 | 1,590 / 1,278 | 3,378 / 1,278 |
| 遺伝子とのリンク（＝発現量の値の数） | 17,870 | 87,130 | 33,717 |
| 用語数 | 37 | 154 | 53 |
| CALOHA / UBERON・CL 対応あり | 37 / 37 | 71 / 101（どちらか 107） | – / 37 |

細胞型クラスは 15。

## 参考

- [sh-ikeda/hpa-expr-spec-rdf-converter](https://github.com/sh-ikeda/hpa-expr-spec-rdf-converter) —
  `proteinatlas.tsv` を RDF 化する gawk スクリプト。組織・細胞型の ID と CALOHA / UBERON / CL 対応表
  （`input/hpa_tissue_caloha.tsv`, `input/hpa_cell_cl.tsv`, `input/hpa_tissue_uberon.tsv`）
