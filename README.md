# ScBRCStates

ScBRCStates applies the breast cancer malignant-cell classifier described in Rott et al. to new scRNA-seq datasets.

The input is an `AnnData` object containing malignant cells. Genes are matched using `adata.var_names`. Model genes absent from the dataset are assigned an expression value of 0 before the same rank transformation used during classifier development.

## Installation

Clone the repository and install it from the repository root:

```bash
pip install .
```

For development/testing:

```bash
pip install -e ".[test]"
```

## Add the classifier

Copy the classifier used in the paper to:

```text
src/scbrcstates/models/scbrcstates_hgb.joblib
```

The joblib file is expected to contain the same objects used in the analysis:

```python
{
    "model": ...,
    "label_encoder": ...,
    "genes": ...,
    "config": ...,
}
```

## Python usage

```python
import scanpy as sc
from scbrcstates import annotate

adata = sc.read_h5ad("malignant_cells.h5ad")
annotate(adata)

print(adata.obs["predicted_cluster"].value_counts())
adata.write_h5ad("malignant_cells_annotated.h5ad")
```

By default the expression matrix is read from `adata.X`. A layer can also be specified:

```python
annotate(adata, layer="log1p")
```

A different output column can be used with:

```python
annotate(adata, key_added="ScBRCStates")
```

## Command line

The same annotation can be run directly on an `.h5ad` file:

```bash
scbrcstates malignant_cells.h5ad malignant_cells_annotated.h5ad
```

Optional arguments:

```bash
scbrcstates input.h5ad output.h5ad --key ScBRCStates
scbrcstates input.h5ad output.h5ad --layer log1p
```

## Input

ScBRCStates is intended for malignant breast cancer cells. It does not identify malignant cells itself.

`adata.var_names` must contain unique gene names. Genes that are not part of the classifier are ignored. Classifier genes absent from the input dataset are kept at 0 before ranking.

The prediction preprocessing is intentionally kept identical to the original analysis:

1. keep model genes present in the dataset;
2. add absent model genes with expression 0;
3. restore the original model-gene order;
4. calculate cell-wise ranks with `np.argsort(np.argsort(...))`;
5. normalize ranks by the maximum rank in each cell;
6. apply the stored HGB classifier and label encoder.

## Tests

Run:

```bash
pytest
```

The tests check that the packaged preprocessing reproduces the original prediction code, including reordered genes and missing genes, and that predictions are correctly added to `adata.obs`.
