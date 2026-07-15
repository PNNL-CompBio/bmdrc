import pandas as pd
import numpy as np
import pytest
import bmdrc
from bmdrc import ContinuousClass


## How to calculate coverage (from within main package directory): 
# coverage run --source=bmdrc -m pytest -x tests/* 

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)

# Continuous Data test
def _continuous():
    obj = ContinuousClass.ContinuousClass(
        df=pd.read_csv("data/Continuous.txt", sep="\t"),
        chemical="Chemical ID",
        concentration="Concentration_uM",
        endpoint="Endpoint",
        response="Measurement"
    )
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    return obj

# Test continous model fitting with a filtered endpoint present
def test_continuous_model_fitting(capsys):
    obj = _continuous()
    obj.fit_models()
    assert not obj.bmds.empty
    assert obj.bmds_filtered is not None

# Test continuous fitting in diagnostic mode
def test_continuous_diagnostic_mode(capsys):
    obj = ContinuousClass.ContinuousClass(
        df=pd.read_csv("data/Continuous.txt", sep="\t"),
        chemical="Chemical ID",
        concentration="Concentration_uM",
        endpoint="Endpoint",
        response="Measurement"
    )
    obj.filter_correlation_score(apply=True)
    obj.fit_models(diagnostic_mode=True)
    assert not obj.bmds.empty

# Test response curve for continuous data
# Runs real fitted curve for the first endpoint, then checks invalid-input validation
def test_continuous_response_curve():
    obj = _continuous()
    ep = obj.bmds["bmdrc.Endpoint.ID"].iloc[0]
    model = str(obj.bmds["Model"].iloc[0]).lower()
    chem = ep.split(" ")[0]
    endpoint= ep.split(" ")[-1]
    obj.response_curve(chemical_name=chem, endpoint_name=endpoint, model=model,
                       fixed_intercept=0)
    # Run curve for every valid continuous model so each model class's 
    # curve-building / predict_x branch is exercised.

    for m in ["power", "hill", "gompertz", "michaelis-mentin", "asymptotic", "weibull"]:
        obj.response_curve(chemical_name=chem, endpoint_name=endpoint, model=m, fixed_intercept=0)
    with pytest.raises(ValueError):
        obj.response_curve(chemical_name="nope", endpoint_name=endpoint, model=model, fixed_intercept=0)
    with pytest.raises(ValueError):
        obj.response_curve(chemical_name=chem, endpoint_name="nope", model=model, fixed_intercept=0)
    with pytest.raises(ValueError):
        obj.response_curve(chemical_name=chem, endpoint_name=endpoint, model="not_a_model", fixed_intercept=0)


# Test continuous benchmark dose and dose table outputs
def test_continuous_outputs(tmp_path):
    obj = _continuous()
    obj.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    obj.output_dose_table(path=str(tmp_path / "dose.csv"))
    assert not pd.read_csv(str(tmp_path / "bmd.csv")).empty

    
#Internal re-fit fails to converge error
def test_continuous_fits_table(tmp_path):
    obj = _continuous()
    obj.output_fits_table(fixed_intercept=0, path=str(tmp_path / "file.csv"))
    assert obj.output_fits_table is not None

