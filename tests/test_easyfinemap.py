"""Tests for EasyFinemap."""

import os
import shutil
from pathlib import Path
import pytest

import pandas as pd

from easyfinemap.loci import Loci
from easyfinemap.constant import ColName
from easyfinemap.easyfinemap import EasyFinemap

PWD = os.path.dirname(os.path.abspath(__file__))
CWD = os.getcwd()


class TestEasyFinemap:
    """Tests for the EasyFinemap class."""

    def test_init(self):
        """Test the EasyFinemap class."""
        if os.path.exists(f"{CWD}/tmp/easyfinemap"):
            shutil.rmtree(f"{CWD}/tmp/easyfinemap")
        easyfinemap = EasyFinemap()
        assert str(easyfinemap.tmp_root) == f"{CWD}/tmp/easyfinemap"
        for file in Path(f"{PWD}/exampledata/").glob("*loci.txt"):
            if file.is_file():
                os.remove(file)
        for file in Path(f"{PWD}/exampledata/").glob("*leadsnp.txt"):
            if file.is_file():
                os.remove(file)


def test_susie_null_n_and_native_cs_merge():
    """SuSiE supports n=NULL and preserves effect-specific CS membership."""
    assert EasyFinemap._format_susie_n_clause(None) == ""
    assert EasyFinemap._format_susie_n_clause(243552.81) == "n=243553,"

    finemap_res = pd.DataFrame(
        {
            ColName.SNPID: ["1-100-A-G", "1-200-C-T"],
            ColName.PP_SUSIE: [0.8, 0.4],
        }
    )
    cs_members = pd.DataFrame(
        {
            ColName.SNPID: ["1-100-A-G", "1-100-A-G", "1-200-C-T"],
            "SUSIE_EFFECT": [1, 2, 2],
            "SUSIE_CS": ["L1", "L2", "L2"],
            "SUSIE_ALPHA": [0.8, 0.3, 0.6],
            "SUSIE_CUM_ALPHA": [0.8, 0.3, 0.9],
        }
    )
    cs_summary = pd.DataFrame(
        {
            "SUSIE_EFFECT": [1, 2],
            "SUSIE_CS": ["L1", "L2"],
            "SUSIE_COVERAGE": [0.95, 0.95],
            "SUSIE_CS_SIZE": [1, 2],
            "SUSIE_MIN_ABS_CORR": [1.0, 0.8],
            "SUSIE_MEAN_ABS_CORR": [1.0, 0.9],
            "SUSIE_MEDIAN_ABS_CORR": [1.0, 0.9],
        }
    )

    merged = EasyFinemap._merge_susie_native_cs(
        finemap_res,
        cs_members,
        cs_summary,
    )

    assert merged.shape[0] == 3
    assert set(merged.loc[merged[ColName.SNPID] == "1-100-A-G", "SUSIE_CS"]) == {
        "L1",
        "L2",
    }
