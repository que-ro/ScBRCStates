from pathlib import Path
import joblib
import numpy as np
from scipy import sparse


MODEL_NAME = "scbrcstates_hgb.joblib"


def _default_model_path():
    return Path(__file__).parent / "models" / MODEL_NAME


def _load_artifacts(model_path=None):
    path = Path(model_path) if model_path else _default_model_path()

    if not path.is_file():
        raise FileNotFoundError(
            f"Classifier not found: {path}\n"
            f"Copy the publication model to src/scbrcstates/models/{MODEL_NAME}."
        )

    artifacts = joblib.load(path)
    required = {"model", "label_encoder", "genes"}
    missing = required.difference(artifacts)
    if missing:
        raise KeyError(f"Missing key(s) in classifier file: {', '.join(sorted(missing))}")

    return artifacts


def _prepare_batch(adata, model_genes, start, stop, layer=None):
    var_names = [str(gene) for gene in adata.var_names]
    if len(var_names) != len(set(var_names)):
        raise ValueError("adata.var_names must contain unique gene names.")

    gene_index = {gene: i for i, gene in enumerate(var_names)}
    model_genes = [str(gene) for gene in model_genes]

    # Start from zeros so genes absent from the dataset keep expression 0.
    X = np.zeros((stop - start, len(model_genes)), dtype=np.float64)

    model_pos = []
    adata_pos = []
    for i, gene in enumerate(model_genes):
        j = gene_index.get(gene)
        if j is not None:
            model_pos.append(i)
            adata_pos.append(j)

    source = adata.X if layer is None else adata.layers[layer]
    if adata_pos:
        block = source[start:stop, adata_pos]
        if sparse.issparse(block):
            block = block.toarray()
        X[:, model_pos] = np.asarray(block, dtype=np.float64)

    # Same ranking used for training and in the original prediction script.
    X = np.argsort(np.argsort(X, axis=1), axis=1)
    X = X / X.max(axis=1, keepdims=True)
    return X


def annotate(
    adata,
    model_path=None,
    key_added="predicted_cluster",
    layer=None,
    batch_size=10000,
    verbose=True,
):
    """Annotate malignant breast cancer cells with ScBRCStates."""
    artifacts = _load_artifacts(model_path)
    model = artifacts["model"]
    label_encoder = artifacts["label_encoder"]
    model_genes = list(artifacts["genes"])

    if batch_size <= 0:
        raise ValueError("batch_size must be greater than 0.")

    available = set(map(str, adata.var_names))
    n_found = sum(gene in available for gene in model_genes)

    if verbose:
        print(f"ScBRCStates: {adata.n_obs:,} cells")
        print(f"Model genes found: {n_found}/{len(model_genes)}")

    predictions = []
    for start in range(0, adata.n_obs, batch_size):
        stop = min(start + batch_size, adata.n_obs)
        X = _prepare_batch(adata, model_genes, start, stop, layer=layer)
        y = model.predict(X)
        predictions.append(label_encoder.inverse_transform(y))

    adata.obs[key_added] = np.concatenate(predictions)

    if verbose:
        print(f"Predictions stored in adata.obs['{key_added}']")

    return adata
