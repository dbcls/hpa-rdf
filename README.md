# HPA RNA specificity → RDF (hpao)

Human Protein Atlas の `proteinatlas.tsv`（最新リリースのみ対応、現在 v25.1）から、
RNA tissue / single cell type specificity を RDF 化する。

## ファイル

| パス | 内容 |
|---|---|
| `hpa2rdf.py` | 変換スクリプト（Python 3 標準ライブラリのみ） |
| `hpao.ttl` | オントロジー（クラス・プロパティ・カテゴリ概念） |
| `mapping/tissue_caloha.tsv`, `mapping/cell_type_caloha.tsv` | HPA 名 → CALOHA の対応（`match` 列: exact / close） |
| `config/hpa/` | rdf-config（`model.yaml`, `prefix.yaml`, `schema.svg`） |
| `sparqlist/*.md` | TogoDX 用 SPARQList 案（`HPA_GRAPH` を読み込み先グラフに置換） |
| `data/proteinatlas.tsv.zip` | 入力（ダウンロードしたもの） |
| `out/hpa.ttl` | 出力 |

## 使い方

```sh
curl -L -o data/proteinatlas.tsv.zip https://www.proteinatlas.org/download/proteinatlas.tsv.zip
python3 hpa2rdf.py data/proteinatlas.tsv.zip -o out/hpa.ttl --hpa-version 25.1
# out/hpa.ttl と hpao.ttl をトリプルストアに読み込む
```

新しいリリースで分類値や列名が変わった場合はエラーで止まる。
対応表にない組織・細胞型名は警告を出し、CALOHA 対応なしの用語として出力する。

## データモデル

```turtle
@prefix hpao: <http://med2rdf.org/hpao/ontology#> .
@prefix hpat: <http://med2rdf.org/hpao/term/> .
@prefix hpa:  <http://med2rdf.org/hpao/resource/> .

ensembl:ENSG00000000003 a obo:SO_0000704 ;
    hpao:tissueSpecificity hpao:LowTissueSpecificity ;
    hpao:tissueDistribution hpao:DetectedInMany ;
    hpao:cellTypeSpecificity hpao:CellTypeGroupEnriched ;
    hpao:cellTypeDistribution hpao:DetectedInMany ;
    hpao:cellTypeSpecificityScore "5"^^xsd:decimal ;
    hpao:groupEnrichedIn hpat:celltype_Esophageal_apical_cells ,
        hpat:celltype_Late_spermatids ;
    hpao:rnaExpression [ hpao:anatomicalEntity hpat:celltype_Esophageal_apical_cells ; hpao:nCPM 1280.5 ] ,
        [ hpao:anatomicalEntity hpat:celltype_Late_spermatids ; hpao:nCPM 3513.3 ] .

hpat:celltype_Late_spermatids a hpao:SingleCellType ;
    rdfs:label "Late spermatids" ;
    skos:closeMatch caloha:TS-2535 .
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
| 用語数（うち CALOHA 対応あり） | 37 (37) | 154 (71) |

## 経緯

TogoDX/Human の旧グラフ（`hpa_tissue_specificity` / `hpa_cell_specificity`、HPA v20 から作成、
`refexo:isPositivelySpecificTo` のみ）を逆算した結果をもとに、分類・分布・スコアを保持するスキーマへ置き換えた。
`mapping/` の CALOHA 対応のうち v20 にあった名前は旧グラフから復元したもの。
