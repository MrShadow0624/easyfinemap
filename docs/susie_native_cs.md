# SuSiE-RSS: official algorithm, signed LD and native credible sets

This implementation delegates inference to **official susieR 0.12.35**. Version-matched algorithm source and package documentation are the primary evidence; HERMES and PolyFun are workflow comparisons. Validation scripts and results are in the parent project's `03_finemap/code/validation/` and `03_finemap/validation/`.

## Official evidence and optional N

The pinned official source is [susie_rss.R at da24d5f](https://github.com/stephenslab/susieR/blob/da24d5fba95845ad082b893f9c80b85b7523c6f6/R/susie_rss.R). Installed `susie_rss`, `susie_get_cs` and `susie_get_pip` formals and function bodies were compared against the downloaded fixed source, not just the package version string.

In this release `missing(n)` selects the no-N branch. The Python API therefore omits the R `n` argument when `sample_size=None`; literal `n=NULL` fails in this version. Providing N selects the finite-sample branch and its PVE-adjusted Z statistics. No-N assumes large N and small effects; the official code recommends providing N when available. These two modes are supported, but are not statistically interchangeable. See the [official summary-statistics tutorial](https://stephenslab.github.io/susieR/articles/finemapping_summary_statistics.html). Current website documentation may describe newer APIs; the fixed source determines this runtime's behavior.

FINEMAP and conditional analysis continue to require N. A factor equivalent N derived from GenomicSEM SE/MAF is not automatically substituted into the no-N analysis. Its applicability depends on factor scaling and the summary-statistics model.

## Signed LD and alignment

RSS requires a signed correlation matrix in the same variant order and effect-allele direction as Z. Squared LD is not a replacement for R. The [official diagnostic tutorial](https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html) demonstrates the consequences of allele inconsistency.

`LDRef.make_ld` now exports BIM A1 dosage with PLINK `--keep-allele-order --recode A`, imputes missing dosages by the per-variant mean, and computes Pearson correlations. See [PLINK's format documentation](https://www.cog-genomics.org/plink/1.9/data#recode). Mean imputation preserves a common sample matrix and positive semidefiniteness. All-missing or zero-variance variants fail explicitly. `intersect` reorders GWAS rows to the extracted BIM and checks chromosome, position and allele pair. `prepare_ld_matrix` applies `D R D` to match GWAS EA; beta and Z retain their original directions.

The generated `.ld.snps.tsv` records SNP order and alleles. Direct callers of `run_susie` must supply an already aligned signed LD matrix; if the sidecar exists it is checked. An unlabelled matrix's biological order/direction cannot be inferred from its dimensions. Forward-strand provenance, genome build and palindromic ambiguity remain input QC responsibilities. Shared LD handling also affects other LD-based methods; this validation targets SuSiE, not all supported tools.

## Native per-effect CSs and achieved coverage

[Fixed-version susie_utils.R](https://github.com/stephenslab/susieR/blob/da24d5fba95845ad082b893f9c80b85b7523c6f6/R/susie_utils.R) defines `susie_get_cs` and `susie_get_pip`. The former uses each effect's alpha row, filters inactive/duplicate/impure sets and returns effect identity in `cs_index`. Marginal PIP is a separate quantity and does not define a multi-effect CS by cumulative locus PIP.

The wrapper uses official membership and purity calculations, with coverage `.95` and minimum absolute correlation `.5` by default. It sets `n_purity` to the number of variants so purity is calculated on each complete CS. The official default samples at most 100 variants in large CSs. L retains the API/CLI default 1; project commands explicitly request L=5. The wrapper uses project `max_iter=200` and official `tol=1e-3`.

A reproducible 0.12.35 metadata bug appears when an earlier CS is removed by purity filtering: returned `coverage` can refer to the wrong effect. For alpha rows `(0.5,0.5,0)` and `(0.96,0.02,0.02)` with identity LD, only effect 2 survives, but the official metadata reports 1 instead of .96. The export computes achieved coverage as `sum(alpha[effect, members])`; it leaves the official inference, membership and purity unchanged. The saved fit's `sets$coverage` is corrected consistently.

CS membership is long format: one SNP can occur in several distinct retained sets. Columns are `SUSIE_EFFECT`, `SUSIE_CS`, `SUSIE_ALPHA`, `SUSIE_CUM_ALPHA`; set summaries add achieved `SUSIE_COVERAGE`, size and min/mean/median absolute correlation. This differs from a single-membership field that silently loses overlaps.

## Diagnostics and retained results

With `output_prefix`, the API retains inputs, signed LD/sidecar, `.fit.rds`, `.variant_pip.tsv`, `.cs_members.tsv`, `.cs_summary.tsv`, `.z_ld_diagnostic.tsv` and `.locus_summary.tsv`. CLI output automatically supplies a per-locus/method prefix under `<outfile>.loci/`. The locus summary records N mode, convergence/iterations, parameters, susieR version, Z–LD consistency and CS counts. A zero-CS locus retains all-variant PIP and stable empty CS tables.

The diagnostic uses official `estimate_s_rss` and `kriging_rss`; `s>.2` or `logLR>2` with `abs(Z)>2` flags project review. The `.2` cutoff is a project convention, not a universal calibrated threshold. It does not automatically flip alleles or modify LD. With external reference LD, `estimate_residual_variance=FALSE`. A nonconverged fit is saved with failure status and raises an error; old successful result files are removed when rerunning a prefix. No successful PIP/CS output is written for that failed fit.

The original arbitrary top-5000-by-P truncation and replacement of infinite effects by +/-100 were removed. Upstream generic CS logic remains for other methods; they have not been validated here.

## HERMES / PolyFun comparison: secondary evidence

The [HERMES wrapper at 486fbcf](https://github.com/ihi-comp-med/hermes2-gwas/blob/486fbcf9ace1c09135bce982f6e6ad7ae67178b7/workflow/rules/scripts/finemap/run_polyfun_susie.py) passes `--method susie` **and `--n`** (floor of maximum locus N). It supports a PolyFun+SuSiE workflow, but provides no evidence for HERMES using no-N.

The [HERMES aggregation script](https://github.com/ihi-comp-med/hermes2-gwas/blob/486fbcf9ace1c09135bce982f6e6ad7ae67178b7/workflow/rules/scripts/finemap/aggregate_fm_locus.R) retains native `CREDIBLE_SET>0`; its cumulative PIP is descriptive rather than a membership criterion.

The [HERMES README](https://github.com/ihi-comp-med/hermes2-gwas/blob/486fbcf9ace1c09135bce982f6e6ad7ae67178b7/README.md) lists PolyFun `v2023-11-14`. [PolyFun finemapper.py at 569f852](https://github.com/omerwe/polyfun/blob/569f852680c9fa707dcdccb3092f0d1f90f1455a/finemapper.py), the last change to that file before that date, obtains PIP with `susie_get_pip` and memberships from native sets. Its single `CREDIBLE_SET` field keeps the first set when SNPs overlap. This date-matched source comparison is not proof of HERMES's exact deployed commit; no exact lock was established. The earlier reference to 7227ed5 (2020) was not sufficient to identify HERMES's runtime version.

## Validation scope

Tests use real rpy2, official R functions in independent processes and PLINK, including no-N/supplied-N, nonuniform prior weights, signed-R allele-flip invariance, native overlap, filtered-effect coverage, zero CSs, nonconvergence and input errors. Full CLI trials use actual ResidualHF_CAD locus data, with exclusions recorded. Nonuniform synthetic weights verify API transmission; they are not actual PolyFun annotations. Full functional-prior generation, all 21 project loci, window/L sensitivity and empirical coverage calibration are separate analyses.
