import anndata as ad

from scbrcstates import annotate


def test_classifier_runs():
    adata = ad.read_h5ad("tests/data/example_adata.h5ad")
    annotate(adata)
    assert "predicted_cluster" in adata.obs.columns
    assert adata.obs["predicted_cluster"].notna().all()
    assert len(adata.obs["predicted_cluster"]) == adata.n_obs


def test_classifier_labels():
    adata = ad.read_h5ad("tests/data/example_adata.h5ad")
    annotate(adata)
    expected_labels = {
        "0",
        "1_4_6",
        "2",
        "3",
        "5",
    }
    observed_labels = set(adata.obs["predicted_cluster"])
    assert observed_labels.issubset(expected_labels)