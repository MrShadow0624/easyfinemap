# SuSiE-RSS: optional sample size and native credible sets

This fork keeps the upstream EasyFinemap workflow, but changes the SuSiE-specific handling in two places.

## 1. Optional sample size in susie_rss

The official SuSiE-RSS method supports analysis without a supplied sample size.
This fork's `environment.yml` pins susieR 0.12.35. In that release, the
no-N branch is activated when the `n` argument is **omitted**; passing the
literal R object `NULL` is not the same as a missing argument in that older
implementation. Therefore, when `--sample-size` is omitted, this fork calls
`susie_rss()` without an `n` argument. When a sample size is supplied, it
is passed normally.

This is the version-compatible implementation of the requested "N = null /
no-N" choice. It preserves the official large-sample/small-effect no-N mode
while remaining compatible with EasyFinemap's pinned susieR version.

Official references:

- current susieR documentation: https://stephenslab.github.io/susieR/reference/susie_rss.html
- susieR 0.12.35 source used by EasyFinemap's environment (no-N branch when `n` is missing): https://github.com/stephenslab/susieR/blob/da24d5fba95845ad082b893f9c80b85b7523c6f6/R/susie_rss.R
- Wang G, Sarkar A, Carbonetto P, Stephens M. A simple new approach to variable selection in regression, with application to genetic fine-mapping. JRSS-B (2020). https://doi.org/10.1111/rssb.12388

## 2. Native, effect-specific SuSiE credible sets

SuSiE is a Sum of Single Effects model. For each effect l, `res$alpha[l, ]` contains effect-specific posterior inclusion probabilities. The official `susie_get_cs()` implementation constructs a credible set separately for each effect and can filter sets by within-set LD purity.

This fork calls `susie_get_cs(res, Xcorr=ld, coverage=..., min_abs_corr=...)` and exports:

- `SUSIE_EFFECT`
- `SUSIE_CS`
- `SUSIE_ALPHA`
- `SUSIE_CUM_ALPHA`
- `SUSIE_COVERAGE`
- `SUSIE_CS_SIZE`
- `SUSIE_MIN_ABS_CORR`
- `SUSIE_MEAN_ABS_CORR`
- `SUSIE_MEDIAN_ABS_CORR`

A variant is allowed to appear in more than one effect-specific credible set. Long-format output preserves those memberships.

Official SuSiE credible-set source:

- https://github.com/stephenslab/susieR/blob/master/R/susie_get_functions.R

## 3. Difference from upstream EasyFinemap

Upstream EasyFinemap v0.4.7 obtains `res$pip` from SuSiE, then its generic `get_credset()` ranks marginal posterior probabilities and cumulatively sums them after multiplying the requested threshold by `max_causal`.

That generic cumulative-marginal-PIP construction remains unchanged for other fine-mapping methods, but is no longer used for `susie` or `polyfun_susie` in this fork. For SuSiE, native per-effect credible sets are used.

Upstream source:

- https://github.com/Jianhua-Wang/easyfinemap/blob/955ed4b4f67ccc469af5a4be3c05fae9002141be/easyfinemap/easyfinemap.py

## 4. HERMES precedent

The HERMES2 heart-failure GWAS used PolyFun + SuSiE for fine-mapping.

Its public wrapper invokes PolyFun with `--method susie`:

- https://github.com/ihi-comp-med/hermes2-gwas/blob/486fbcf9ace1c09135bce982f6e6ad7ae67178b7/workflow/rules/scripts/finemap/run_polyfun_susie.py

The PolyFun implementation used in that workflow obtains marginal PIP with `susie_get_pip(susie_obj)`, but assigns `CREDIBLE_SET` from native `susie_obj$sets` rather than reconstructing sets from marginal PIP:

- https://github.com/omerwe/polyfun/blob/7227ed5261ee0ed9ad1af1031419a2180099b953/finemapper.py

The HERMES aggregation script retains variants with `CREDIBLE_SET > 0`. Its `cum_pip` column is descriptive and does not define credible-set membership:

- https://github.com/ihi-comp-med/hermes2-gwas/blob/486fbcf9ace1c09135bce982f6e6ad7ae67178b7/workflow/rules/scripts/finemap/aggregate_fm_locus.R

Paper:

- Shah S et al. Genome-wide association study meta-analysis provides insights into the etiology of heart failure and its subtypes. Nature Genetics. https://doi.org/10.1038/s41588-024-02064-3

## Scope

These changes intentionally do not alter FINEMAP, PAINTOR, CAVIARBF, ABF, locus definition, conditional analysis, or PolyFun prior generation. The goal is only to make SuSiE-RSS sample-size handling explicit and to report SuSiE credible sets according to the SuSiE model itself.
