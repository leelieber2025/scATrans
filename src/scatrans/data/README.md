# Data Folder

This folder contains gene feature files used for bias correction.

## Files

- `mouse_2020A_gene_features.parquet`  
  Precomputed gene features (gene_length + intron_number) for mouse genes.  
  Used automatically by `active_score()` for Huber bias correction when present in `adata.var`.

You can generate similar files for other species using:
```python
from scatrans.pp_bias import generate_gene_features_from_gtf
```
