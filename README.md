# easyfinemap


[![pypi](https://img.shields.io/pypi/v/easyfinemap.svg)](https://pypi.org/project/easyfinemap/)
[![python](https://img.shields.io/pypi/pyversions/easyfinemap.svg)](https://pypi.org/project/easyfinemap/)
[![Build Status](https://github.com/Jianhua-Wang/easyfinemap/actions/workflows/dev.yml/badge.svg)](https://github.com/Jianhua-Wang/easyfinemap/actions/workflows/dev.yml)
[![codecov](https://codecov.io/gh/Jianhua-Wang/easyfinemap/branch/main/graphs/badge.svg)](https://codecov.io/github/Jianhua-Wang/easyfinemap)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI download month](https://img.shields.io/pypi/dm/easyfinemap.svg)](https://pypi.org/project/easyfinemap/)
[![Build Status](https://github.com/Jianhua-Wang/easyfinemap/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/Jianhua-Wang/easyfinemap/actions/workflows/python-package-conda.yml)


A flexible framework for GWAS fine-mapping


* Documentation: <https://Jianhua-Wang.github.io/easyfinemap>
* GitHub: <https://github.com/Jianhua-Wang/easyfinemap>
* PyPI: <https://pypi.org/project/easyfinemap/>
* License: MIT


## Features
* Formatting of summary statistics using smunger.
* Fast extraction of summary statistics for fine-mapping using tabix/bgzip.
* Checking and formatting of LD reference.
* Representation of variants using Unique SNP ID.
* Support for identifying independent loci using three methods: distance, LD clumping, and conditional analysis.
* Fine-mapping without the need for LD reference.
* Support for four LD-based fine-mapping tools and function-based fine-mapping.
* Fine-mapping combined with conditional analysis.

<!-- ## Finemapping approaches
* LD-free
    * aBF
* LD-based
    * FINEMAP
    * CAVIARBF
    * PAINTOR -->

## Fork-specific SuSiE-RSS behavior

This fork uses official susieR with explicit inputs and native CS outputs:

1. `--sample-size` is optional for SuSiE-RSS. If omitted, the underlying
   `susieR::susie_rss()` call omits the `n` argument (the official no-N
   path in the fork's pinned susieR 0.12.35); if supplied, the numerical
   sample size is passed through.
2. SuSiE credible sets are taken from native, effect-specific
   `susie_get_cs()` output (including coverage and LD-purity summaries),
   rather than being reconstructed by cumulatively summing marginal PIPs
   across a locus.
3. LD uses signed dosage correlations, reordered to BIM and aligned to GWAS
   effect alleles. Squared LD (`r²`) is not passed as the RSS correlation matrix.
4. CLI runs retain each fit, all-variant PIP, native CS tables, diagnostics and
   locus status under `<outfile>.loci/`, including loci without retained CSs.
   Nonconverged fits raise an error after saving a checkpoint.

For the project preset, pass `--max-causal 5 --credible-threshold 0.95
--susie-min-abs-corr 0.5` explicitly; the CLI's L default remains 1.

Implementation notes and references are documented in
[`docs/susie_native_cs.md`](docs/susie_native_cs.md).
