"""用真实 susieR/rpy2/PLINK 验证算法结果与输入对齐，不替换算法为 mock。
在已安装 susieR 0.12.35 的环境运行；合成输入和结果写入 pytest 临时目录。
"""
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import rpy2.robjects as ro
from rpy2.rinterface_lib.embedded import RRuntimeError

from easyfinemap.easyfinemap import EasyFinemap


@pytest.mark.parametrize("sample_size", [0, 1, np.nan, np.inf])
def test_invalid_n(sample_size):
    with pytest.raises(ValueError, match="sample_size"):
        EasyFinemap._format_susie_n_clause(sample_size)


def test_missing_n_differs_from_null_in_pinned_version():
    if str(ro.r("as.character(packageVersion('susieR'))")[0]) == "0.12.35":
        with pytest.raises(RRuntimeError):
            ro.r("susieR::susie_rss(z=c(4,3),R=diag(2),n=NULL,L=1)")


def test_signed_r_and_allele_flip_invariance(tmp_path):
    matrix = np.array([[1., -.8], [-.8, 1.]])
    d = pd.DataFrame({"SNPID": ["a", "b"], "BETA": [8., -6.4], "SE": [1., 1.]})
    ld = tmp_path / "ld"
    np.savetxt(ld, matrix)
    correct = EasyFinemap().run_susie(d, str(ld), max_causal=2)
    flipped = d.copy()
    flipped.loc[1, "BETA"] *= -1
    np.savetxt(ld, matrix * np.outer([1, -1], [1, -1]))
    consistent = EasyFinemap().run_susie(flipped, str(ld), max_causal=2)
    np.testing.assert_allclose(correct, consistent, atol=1e-12, rtol=0)
    np.savetxt(ld, matrix ** 2)
    incorrect = EasyFinemap().run_susie(d, str(ld), max_causal=2)
    assert correct.iloc[1] < .001 and incorrect.iloc[1] > .99


@pytest.mark.parametrize("sample_size", [None, 1000])
@pytest.mark.parametrize("weighted", [False, True])
def test_matches_official_r(tmp_path, sample_size, weighted):
    # Step1. 构造含负相关的 LD；非均匀先验只测试权重传递，不代表真实 PolyFun 先验。
    matrix = np.array([[1., -.8, 0.], [-.8, 1., 0.], [0., 0., 1.]])
    ld_file = tmp_path / "signed.ld"
    np.savetxt(ld_file, matrix)
    d = pd.DataFrame({"SNPID": ["1-100-A-G", "1-200-C-T", "1-300-A-C"], "EA": ["A", "C", "A"],
                      "BETA": [5., -4., 4.5], "SE": [1., 1., 1.], "SNPVAR": [1., 3., 2.]})
    prefix = str(tmp_path / "fit")
    susie_input = d if weighted else d.drop(columns='SNPVAR')
    pip, members, summary = EasyFinemap().run_susie(susie_input, str(ld_file), sample_size=sample_size, max_causal=2,
                         return_native_cs=True, output_prefix=prefix)
    expected_weights = d.SNPVAR / d.SNPVAR.sum() if weighted else np.repeat(1 / len(d), len(d))
    np.testing.assert_allclose(pd.read_csv(f"{prefix}.inputs.tsv", sep="\t").SNPVAR, expected_weights, atol=1e-15, rtol=0)

    # Step2. 在独立 R 进程直接调用官方函数，比较全部 PIP 与每个 CS 的成员和 alpha。
    n_clause = "" if sample_size is None else f"n={sample_size},"
    oracle = tmp_path / "oracle.R"
    oracle.write_text(f"""
library(data.table)
library(susieR)
args = commandArgs(TRUE)
d = fread(paste0(args[1], '.inputs.tsv'))
R = as.matrix(fread(args[2], header=FALSE))
fit = susie_rss(z=d$BETA/d$SE, R=R, {n_clause} L=2, prior_weights=d$SNPVAR,
                estimate_residual_variance=FALSE, coverage=.95, min_abs_corr=.5, n_purity=nrow(d), max_iter=200, tol=1e-3)
saveRDS(fit, paste0(args[1], '.oracle.rds'))
fwrite(data.table(PIP=susie_get_pip(fit)), paste0(args[1], '.oracle.tsv'), sep='\\t')
""")
    subprocess.run(["Rscript", str(oracle), prefix, str(ld_file)], check=True, capture_output=True, text=True)
    np.testing.assert_allclose(pip, pd.read_csv(f"{prefix}.oracle.tsv", sep="\t")["PIP"], rtol=0, atol=1e-12)
    ro.globalenv["oracle_file"] = ro.StrVector([f"{prefix}.oracle.rds"])
    ro.r("oracle_fit = readRDS(oracle_file)")
    assert bool(ro.r("oracle_fit$converged")[0])
    cs = ro.r("oracle_fit$sets$cs")
    effect = np.asarray(ro.r("oracle_fit$sets$cs_index"), dtype=int)
    alpha = np.asarray(ro.r("oracle_fit$alpha"))
    assert len(cs) == len(summary)
    for k in range(len(cs)):  # 验证每个单效应，不能用 locus 总 PIP 代替 coverage。
        idx = np.asarray(cs[k], dtype=int) - 1
        actual = members[members.SUSIE_EFFECT == effect[k]]
        assert set(actual.SNPID) == set(d.iloc[idx].SNPID)
        expected_alpha = dict(zip(d.SNPID, alpha[effect[k] - 1]))
        np.testing.assert_allclose(actual.SUSIE_ALPHA, actual.SNPID.map(expected_alpha), atol=1e-12, rtol=0)
        achieved = summary.loc[summary.SUSIE_EFFECT == effect[k], "SUSIE_COVERAGE"].iloc[0]
        assert achieved == pytest.approx(alpha[effect[k] - 1, idx].sum(), abs=1e-12)


def test_no_cs_preserves_outputs(tmp_path):
    ld = tmp_path / "null.ld"
    np.savetxt(ld, np.eye(3))
    d = pd.DataFrame({"SNPID": ["a", "b", "c"], "BETA": [0., 0., 0.], "SE": [1., 1., 1.]})
    pip, members, summary = EasyFinemap().run_susie(d, str(ld), return_native_cs=True, output_prefix=str(tmp_path / "null"))
    assert members.empty and summary.empty
    merged = EasyFinemap._merge_susie_native_cs(d, members, summary)
    assert merged.empty and "SUSIE_ALPHA" in merged and "SUSIE_COVERAGE" in merged
    assert len(pip) == 3 and np.isfinite(pip).all()
    status = pd.read_csv(tmp_path / "null.locus_summary.tsv", sep="\t").iloc[0]
    assert status.converged and status.n_credible_sets == 0 and status.n_variants_in_cs == 0
    assert (tmp_path / "null.fit.rds").exists()
    assert len(pd.read_csv(tmp_path / "null.variant_pip.tsv", sep="\t")) == 3


def test_nonconvergence_saves_checkpoint(tmp_path):
    ld = tmp_path / "ld"
    np.savetxt(ld, np.eye(3))
    d = pd.DataFrame({"SNPID": ["a", "b", "c"], "BETA": [8., 6., .1], "SE": [1., 1., 1.]})
    (tmp_path / "failed.variant_pip.tsv").write_text("stale_success\n")
    (tmp_path / "failed.cs_members.tsv").write_text("stale_success\n")
    with pytest.raises(RRuntimeError, match="did not converge"):
        EasyFinemap().run_susie(d, str(ld), max_causal=2, susie_max_iter=1, output_prefix=str(tmp_path / "failed"))
    assert (tmp_path / "failed.fit.rds").exists()
    assert not pd.read_csv(tmp_path / "failed.locus_summary.tsv", sep="\t").converged.iloc[0]
    assert not (tmp_path / "failed.variant_pip.tsv").exists()
    assert not (tmp_path / "failed.cs_members.tsv").exists()


def test_native_overlap_and_coverage_index():
    # 官方函数直接确认重叠 CS 和 0.12.35 purity 筛选后的 coverage 索引缺陷。
    ro.r("""
library(susieR)
overlap_fit = list(alpha=rbind(c(.49,.49,.02),c(.02,.49,.49)), V=c(1,1))
overlap_cs = susie_get_cs(overlap_fit, Xcorr=matrix(.8,3,3)+diag(.2,3), coverage=.95, n_purity=3)
filtered_fit = list(alpha=rbind(c(.5,.5,0),c(.96,.02,.02)), V=c(1,1))
filtered_cs = susie_get_cs(filtered_fit, Xcorr=diag(3), coverage=.95, min_abs_corr=.5, n_purity=3)
""")
    assert 2 in np.asarray(ro.r("overlap_cs$cs[[1]]")) and 2 in np.asarray(ro.r("overlap_cs$cs[[2]]"))
    assert list(ro.r("filtered_cs$cs_index")) == [2]
    achieved = float(ro.r("sum(filtered_fit$alpha[filtered_cs$cs_index[1],filtered_cs$cs[[1]]])")[0])
    assert achieved == pytest.approx(.96)
    if str(ro.r("as.character(packageVersion('susieR'))")[0]) == "0.12.35":
        assert float(ro.r("filtered_cs$coverage[1]")[0]) == pytest.approx(1.)


def test_filtered_effect_exports_correct_coverage(tmp_path):
    # 真实拟合也会触发索引缺陷：第一个 effect 的多 SNP CS 不纯，第二个 effect 的 singleton 保留。
    ld = tmp_path / "ld"
    np.savetxt(ld, np.eye(4))
    d = pd.DataFrame({"SNPID": ["a", "b", "c", "d"], "BETA": [6., 6., 4., 0.], "SE": [1., 1., 1., 1.]})
    prefix = str(tmp_path / "filtered")
    pip, members, summary = EasyFinemap().run_susie(d, str(ld), max_causal=3, return_native_cs=True, output_prefix=prefix)
    assert members.SUSIE_EFFECT.tolist() == [2] and members.SNPID.tolist() == ["c"]
    assert summary.SUSIE_COVERAGE.iloc[0] == pytest.approx(members.SUSIE_ALPHA.sum(), abs=1e-12)
    assert .98 < summary.SUSIE_COVERAGE.iloc[0] < .995
    ro.globalenv["filtered_saved_file"] = ro.StrVector([f"{prefix}.fit.rds"])
    ro.r("filtered_saved_fit = readRDS(filtered_saved_file); raw_sets = susieR::susie_get_cs(filtered_saved_fit,Xcorr=diag(4),n_purity=4)")
    assert float(ro.r("filtered_saved_fit$sets$coverage[1]")[0]) == pytest.approx(summary.SUSIE_COVERAGE.iloc[0], abs=1e-12)
    if str(ro.r("as.character(packageVersion('susieR'))")[0]) == "0.12.35":
        assert abs(float(ro.r("raw_sets$coverage[1]")[0]) - summary.SUSIE_COVERAGE.iloc[0]) > .005


@pytest.mark.parametrize("bad", ["se", "prior", "order", "allele"])
def test_invalid_input_fails_before_fit(tmp_path, bad):
    ld = tmp_path / "ld"
    np.savetxt(ld, np.eye(2))
    d = pd.DataFrame({"SNPID": ["a", "b"], "EA": ["A", "C"], "BETA": [4., 3.], "SE": [1., 1.], "SNPVAR": [1., 1.]})
    sidecar = d[["SNPID", "EA"]].copy()
    if bad == "se":
        d.loc[0, "SE"] = 0
    if bad == "prior":
        d.loc[0, "SNPVAR"] = -1
    if bad == "order":
        sidecar = sidecar.iloc[::-1]
    if bad == "allele":
        sidecar.loc[0, "EA"] = "G"
    sidecar.to_csv(f"{ld}.snps.tsv", sep="\t", index=False)
    with pytest.raises(ValueError):
        EasyFinemap().run_susie(d, str(ld), prior_file="preannotated", output_prefix=str(tmp_path / "bad"))
    assert not (tmp_path / "bad.fit.rds").exists()


def test_plink_signed_ld_and_cli(tmp_path):
    # Step1. 真实生成 PLINK BED；已知剂量的两列负相关，同时含一个缺失基因型。
    genotype = np.array([[0, 2, 0], [1, 1, 1], [2, 0, 0], [0, 2, 2], [1, 1, 1], [2, 0, 2], [1, 2, np.nan], [0, 1, 1]], dtype=float)
    pairs = [("A", "G"), ("C", "T"), ("A", "C")]
    ids = ["1-100-A-G", "1-200-C-T", "1-300-A-C"]
    source = tmp_path / "panel"
    source.with_suffix(".map").write_text("".join(f"1 {snp} 0 {bp}\n" for snp, bp in zip(ids, [100, 200, 300])))
    with source.with_suffix(".ped").open("w") as out:
        for i, row in enumerate(genotype):  # PED 是按样本顺序逐行写出的测试格式。
            alleles = []
            for dosage, (a1, a2) in zip(row, pairs):
                alleles.extend(["0", "0"] if np.isnan(dosage) else [a1] * int(dosage) + [a2] * (2 - int(dosage)))
            out.write(f"F I{i} 0 0 1 -9 " + " ".join(alleles) + "\n")
    subprocess.run(["plink", "--file", str(source), "--make-bed", "--keep-allele-order", "--out", str(source)], check=True, capture_output=True)
    # Step2. GWAS 行顺序反转且第二个 SNP 的 EA 反转；直接剂量作为独立 LD 参考。
    d = pd.DataFrame({"CHR": [1, 1, 1], "BP": [100, 200, 300], "SNPID": ids, "rsID": ["rs1", "rs2", "rs3"],
                      "EA": ["A", "T", "A"], "NEA": ["G", "C", "C"], "BETA": [8., 6.4, .1], "SE": [1., 1., 1.], "P": [1e-15, 1e-10, .9]})
    d = d.iloc[::-1].reset_index(drop=True)
    aligned = EasyFinemap().prepare_ld_matrix(d, str(source), str(tmp_path / "aligned"), use_ref_EAF=True)
    imputed = np.where(np.isnan(genotype), np.nanmean(genotype, axis=0), genotype)
    imputed[:, 1] = 2 - imputed[:, 1]
    expected = np.corrcoef(imputed, rowvar=False)
    np.testing.assert_allclose(np.loadtxt(tmp_path / "aligned.ld"), expected, rtol=0, atol=1e-14)
    assert aligned.SNPID.tolist() == ids
    assert aligned.EAF.iloc[1] == pytest.approx((2 - genotype[:, 1]).mean() / 2, abs=1e-6)
    assert np.corrcoef(genotype[:, :2], rowvar=False)[0, 1] < 0

    # Step3. 子进程跑完整 CLI no-N 路径，验证多进程及持久输出。
    d.to_csv(tmp_path / "gwas.tsv", sep="\t", index=False)
    pd.DataFrame({"CHR": [1], "START": [50], "END": [350], "LEAD_SNP": [ids[0]]}).to_csv(tmp_path / "loci.tsv", sep="\t", index=False)
    d.iloc[[2]].to_csv(tmp_path / "lead.tsv", sep="\t", index=False)
    outfile = tmp_path / "cli.tsv"
    result = subprocess.run([sys.executable, "-m", "easyfinemap.cli", "fine-mapping", str(tmp_path / "gwas.tsv"),
               str(tmp_path / "loci.tsv"), str(tmp_path / "lead.tsv"), str(outfile), "-m", "susie", "--ldref", str(source),
               "--max-causal", "2", "--credible-threshold", ".95", "--threads", "1"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    output = pd.read_csv(outfile, sep="\t")
    assert "SUSIE_EFFECT" in output and len(output) > 0
    status = pd.read_csv(f"{outfile}.loci/1_50_350.susie.locus_summary.tsv", sep="\t").iloc[0]
    assert status.converged and status.n_mode == "none_large_sample_small_effect"

    # Step4. 同一真实 LD 面板试跑预注释 SNPVAR；不提供 prior-file 或 tabix 索引。
    d['SNPVAR'] = [2., 3., 1.]
    d.to_csv(tmp_path / "gwas.tsv", sep="\t", index=False)
    weighted_outfile = tmp_path / "polyfun.tsv"
    result = subprocess.run([sys.executable, "-m", "easyfinemap.cli", "fine-mapping", str(tmp_path / "gwas.tsv"),
               str(tmp_path / "loci.tsv"), str(tmp_path / "lead.tsv"), str(weighted_outfile), "-m", "polyfun_susie", "--ldref", str(source),
               "--max-causal", "2", "--credible-threshold", ".95", "--threads", "1"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    weighted_input = pd.read_csv(f"{weighted_outfile}.loci/1_50_350.polyfun_susie.inputs.tsv", sep="\t")
    expected = d.set_index('SNPID').SNPVAR / d.SNPVAR.sum()
    np.testing.assert_allclose(weighted_input.SNPVAR, weighted_input.SNPID.map(expected), atol=1e-15, rtol=0)
    assert "SUSIE_EFFECT" in pd.read_csv(weighted_outfile, sep="\t")
