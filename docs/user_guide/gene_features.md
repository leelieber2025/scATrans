# Gene Feature Attachment & CLI

Gene length and intron count are used for optional bias correction inside
`active_score`.

```python
# Use bundled tables
adata = scat.add_gene_features(adata, organism="mouse")  # or "human"

# or provide your own table
adata = scat.add_gene_features(adata, gene_features_path="my_features.parquet")
```

The package includes tables for mouse and human. Use `organism="mouse"`
(default) or `organism="human"` when calling `add_gene_features`. For other
species or custom annotations use the gene feature generator CLI.

## Command-line interface

The only console script is the gene-feature table generator:

```bash
pip install "scatrans[gene_features]"
generate-gene-features --gtf /path/to/genes.gtf --output my_features.parquet --organism human
```

Works with 10x `genes.gtf` or GENCODE GTFs:

```bash
# Mouse
generate-gene-features --gtf /path/to/genes.gtf \
                       --output my_mouse_features.parquet \
                       --organism mouse

# Human (GENCODE or 10x)
generate-gene-features --gtf gencode.v49.primary_assembly.annotation.gtf \
                       --output human_GRCh38_2024A_gene_features.parquet \
                       --organism human
```

Then use it:

```python
import scatrans as scat

adata = scat.add_gene_features(
    adata,
    gene_features_path="human_GRCh38_2024A_gene_features.parquet"
)

# bias correction will now be able to use length + intron_number
adata_res, significant, all_results = scat.active_score(adata, ...)
```

You can also call the generator programmatically:

```python
from scatrans import generate_gene_features_from_gtf

df = generate_gene_features_from_gtf(
    "path/to/genes.gtf",
    output_name="my_custom_features.parquet",
    organism="human"
)
```

See also `scat.list_available_gene_features()` (for bundled tables) and the
full signature of `add_gene_features` in the {doc}`../api/index`.

**Tip**: The generated parquet must contain a `gene_name` column (plus
`gene_length` and `intron_number`). `add_gene_features` does a `reindex` on
your `adata.var_names`.

**GTF notes**: Exon lines need a `transcript_id` attribute (slim GFF3
conversions often drop it — you will get a clear error). Genes without a
usable exon union get **NaN** length (not 0). Huber bias correction in
`active_score` only fits genes with `gene_length > 0`; missing-length genes
use median residual centering when expressed.
