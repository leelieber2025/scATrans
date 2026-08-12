# Security Policy

## Supported Versions

scATrans is under active development. Security fixes are applied to the
latest released version on [PyPI](https://pypi.org/project/scatrans/); we do
not maintain long-term security support for older releases.

| Version | Supported |
| ------- | --------- |
| Latest release | :white_check_mark: |
| Older releases | :x: |

## Reporting a Vulnerability

Please **do not** open a public GitHub issue for security vulnerabilities.

Instead, use GitHub's private vulnerability reporting:
[Report a vulnerability](https://github.com/leelieber2025/scATrans/security/advisories/new)
(repository → **Security** tab → **Report a vulnerability**).

Include a description of the issue, steps to reproduce, and the affected
version. We will acknowledge reports and follow up with a timeline for a fix
once the report is triaged.

## Scope

scATrans is a local single-cell analysis library (no network services, no
telemetry, no remote code execution by design). Relevant concerns include:
unsafe deserialization of untrusted input files, path traversal in file I/O
helpers (e.g. `save_enrichment_report`, gene-feature generation), and
dependency vulnerabilities in the packages it relies on (scanpy, anndata,
PyDESeq2, gseapy, etc.) — for the latter, please also check whether the
issue originates upstream.
