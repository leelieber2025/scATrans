# scATrans

Single-cell Active Transcription Analysis

## New in this version: Dual-Track Design

- `mode="heuristic"` (default): Original fast method
- `mode="advanced"`: Uses scVelo moments + Huber correction (experimental)

See `tl.py` docstring and examples for usage.

## Installation

```bash
pip install -e ".[advanced]"
```
