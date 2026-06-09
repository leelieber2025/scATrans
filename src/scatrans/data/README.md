# Data Folder

This folder contains gene feature files used for bias correction.

## Files

- `mouse_2020A_gene_features.parquet`  
  Precomputed gene features (gene_length + intron_number) for mouse genes.  
  Used automatically by `active_score()` for Huber bias correction when present in `adata.var`.

You can generate similar files for other species using:
```bash
generate-gene-features --gtf /path/to/genes.gtf --output human_gencode_v49_gene_features.parquet --organism human
```

Or programmatically:
```python
from scatrans.pp_bias import generate_gene_features_from_gtf
generate_gene_features_from_gtf("genes.gtf", output_name="human_features.parquet")
```
