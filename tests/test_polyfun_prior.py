"""验证预注释 SNPVAR、严格 QC 与旧 tabix prior-file 入口。
在 EasyFinemap 测试环境运行；模拟外部文件查询，真实 SuSiE 拟合写入 pytest 临时目录。
"""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from easyfinemap.easyfinemap import EasyFinemap


@pytest.fixture
def prior_sumstats():
    return pd.DataFrame({"CHR": [1, 1, 1], "BP": [100, 200, 300],
                         "SNPID": ["1-100-A-G", "1-200-C-T", "1-300-A-C"], "rsID": ["rs1", "rs2", "rs3"],
                         "EA": ["A", "C", "A"], "NEA": ["G", "T", "C"], "BETA": [5., -4., 4.5],
                         "SE": [1., 1., 1.], "P": [1e-6, 1e-4, 1e-5], "SNPVAR": [1., 4., 2.]})


@pytest.fixture
def aligned_ld(monkeypatch):
    # 这里只替换面板提取；SuSiE 仍使用真实 R 算法，CLI 的真实 PLINK 另有回归测试。
    def prepare_ld(self, sumstats, outprefix, **kwargs):
        np.savetxt(f"{outprefix}.ld", [[1., -.8, 0.], [-.8, 1., 0.], [0., 0., 1.]])
        return sumstats.copy()

    monkeypatch.setattr(EasyFinemap, "prepare_ld_matrix", prepare_ld)


def test_preannotated_polyfun_susie(tmp_path, prior_sumstats, aligned_ld):
    prefix = str(tmp_path / "fit")
    result = EasyFinemap().finemap_locus(prior_sumstats, "1", 50, 350, ["polyfun_susie", "susie"],
             "1-100-A-G", max_causal=2, output_prefix=prefix)
    weighted = pd.read_csv(f"{prefix}.polyfun_susie.inputs.tsv", sep="\t")
    uniform = pd.read_csv(f"{prefix}.susie.inputs.tsv", sep="\t")
    np.testing.assert_allclose(weighted.SNPVAR, [1 / 7, 4 / 7, 2 / 7], atol=1e-15, rtol=0)
    np.testing.assert_allclose(uniform.SNPVAR, [1 / 3] * 3, atol=1e-15, rtol=0)
    assert np.isfinite(result[['PP_POLYFUN_SUSIE', 'PP_SUSIE']]).all().all()
    assert prior_sumstats.SNPVAR.tolist() == [1., 4., 2.]  # 输入对象不能被归一化覆盖。


@pytest.mark.parametrize("bad_prior", [0., -1., np.nan, np.inf, -np.inf, pd.NA, None])
def test_preannotated_prior_rejects_invalid(tmp_path, prior_sumstats, bad_prior):
    prior_sumstats.loc[1, 'SNPVAR'] = bad_prior
    with pytest.raises(ValueError, match="SNPVAR"):
        EasyFinemap().run_susie(prior_sumstats, str(tmp_path / "unused.ld"))


def test_polyfun_susie_requires_prior(prior_sumstats):
    with pytest.raises(ValueError, match="SNPVAR or --prior-file"):
        EasyFinemap().finemap_locus(prior_sumstats.drop(columns='SNPVAR'), "1", 50, 350,
                                   ["polyfun_susie"], "1-100-A-G")


@pytest.fixture
def legacy_prior(tmp_path, monkeypatch):
    prior_file = tmp_path / "prior.tsv.gz"
    pd.DataFrame(columns=['CHR', 'BP', 'A1', 'A2', 'snpvar_bin']).to_csv(prior_file, sep="\t", index=False)
    (tmp_path / "prior.tsv.gz.tbi").touch()
    rows = [['1', '100', 'G', 'A', '1'], ['1', '200', 'C', 'T', '4'], ['1', '300', 'A', 'C', '2']]

    def query_prior(region):
        assert region == '1:100-300'  # legacy BP 查询必须包含闭区间首尾。
        return iter(rows)

    monkeypatch.setattr('easyfinemap.easyfinemap.tabix.open', lambda path: SimpleNamespace(querys=query_prior))
    return str(prior_file), rows


def test_legacy_prior_file_compatible(tmp_path, prior_sumstats, aligned_ld, legacy_prior):
    prior_file, _ = legacy_prior
    prefix = str(tmp_path / "legacy")
    prior_sumstats['SNPVAR'] = 99.  # 显式 prior-file 沿用旧入口，并覆盖预注释列。
    result = EasyFinemap().finemap_locus(prior_sumstats, "1", 50, 350, ["polyfun_susie"],
             "1-100-A-G", prior_file=prior_file, max_causal=2, output_prefix=prefix)
    weighted = pd.read_csv(f"{prefix}.polyfun_susie.inputs.tsv", sep="\t")
    np.testing.assert_allclose(weighted.SNPVAR, [1 / 7, 4 / 7, 2 / 7], atol=1e-15, rtol=0)
    assert np.isfinite(result.PP_POLYFUN_SUSIE).all()


@pytest.mark.parametrize("bad_prior", ['missing', '0', '-1', 'nan', 'inf'])
def test_legacy_prior_rejects_invalid(prior_sumstats, legacy_prior, bad_prior):
    prior_file, rows = legacy_prior
    if bad_prior == 'missing':
        rows.pop()
    else:
        rows[1][-1] = bad_prior
    with pytest.raises(ValueError, match="SNPVAR"):
        EasyFinemap().annotate_prior(prior_sumstats, prior_file)
