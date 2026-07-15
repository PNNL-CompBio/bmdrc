import pandas as pd
import numpy as np
import pytest
from bmdrc import ContinuousClass

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)
# Fit continuous data
def _continuous(apply_filter=True):
    obj = ContinuousClass.ContinuousClass(
        df=pd.read_csv("data/Continuous.txt", sep="\t").drop("Notes", axis = 1),
        chemical="Chemical ID",
        concentration="Concentration_uM",
        endpoint="Endpoint",
        response="Measurement"
    )
    if apply_filter:
        obj.filter_correlation_score(apply=True)
    obj.fit_models()
    return obj

# Test benchmark dose output with a filtered endpoint (correlation score filter)
def test_continuous_benchmark_dose_filtered(tmp_path):
    obj = _continuous(apply_filter=True)
    assert obj.bmds_filtered is not None
    obj.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    written = pd.read_csv(str(tmp_path / "bmd.csv"))
    assert not written.empty
    # Filtered  endpoint should be flagged with the correlation score
    assert (written["Modeled_Flag"] == "Fail - correlation score filter").any()

# Test benchmark dose output when nothing was filtered out.
def test_continuous_benchmark_dose_no_filtered(tmp_path):
    obj = ContinuousClass.ContinuousClass(
        df=pd.read_csv("data/Continuous.txt", sep="\t").drop("Notes", axis = 1),
        chemical="Chemical ID",
        concentration="Concentration_uM",
        endpoint="Endpoint",
        response="Measurement"
    )
    obj.filter_correlation_score(score=-1, apply=True)
    obj.fit_models()
    obj.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    assert not pd.read_csv(str(tmp_path / "bmd.csv")).empty

# Test dose table output
def test_continuous_dose_table(tmp_path):
    obj = _continuous(apply_filter=True)
    obj.output_dose_table(path=str(tmp_path / "dose.csv"))
    assert not pd.read_csv(str(tmp_path / "dose.csv")).empty

# Test benchmark dose output with no path
def test_continuous_benchmark_dose_no_path():
    obj = _continuous(apply_filter=True)
    obj.output_benchmark_dose(path=None)

# Test continuous benchmark dose rebuilds plate groups
def test_continuous_benchmark_dose_rebuilds_plate_groups(tmp_path):
    obj = _continuous()
    del obj.plate_groups
    obj.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    assert hasattr (obj, "plate_groups")