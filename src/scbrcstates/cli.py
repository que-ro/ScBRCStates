import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Annotate malignant breast cancer cells with ScBRCStates."
    )
    parser.add_argument("input", help="Input .h5ad file")
    parser.add_argument("output", help="Output .h5ad file")
    parser.add_argument("--key", default="predicted_cluster", help="Column added to adata.obs")
    parser.add_argument("--layer", default=None, help="AnnData layer to use instead of adata.X")
    parser.add_argument("--model", default=None, help="Path to a classifier .joblib file")
    args = parser.parse_args()

    import anndata as ad
    from .classifier import annotate

    adata = ad.read_h5ad(args.input)
    annotate(
        adata,
        model_path=args.model,
        key_added=args.key,
        layer=args.layer,
    )
    adata.write_h5ad(args.output)


if __name__ == "__main__":
    main()
