import pandas as pd
import numpy as np
import pytest
import bmdrc
from bmdrc import BinaryClass
from bmdrc import ProportionalClass
from bmdrc.model_fitting import Calculate_BMD, Calculate_BMDL, Logistic
from bmdrc.model_fitting import Logistic

### Must use pandas version < 3, works with pandas 2.3.3 v. 
## How to calculate coverage (from within main package directory): 
# coverage run --source=bmdrc -m pytest -x tests/* 

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)


# Test binary class model fitting
def _binary():
    obj = BinaryClass.BinaryClass(
        df=pd.read_csv("data/Binary_Simplified_Long.csv"),
        chemical="chemical.id",
        plate="plate.id",
        well="well",
        concentration="concentration",
        endpoint="endpoint",
        value="value",
        format="long"
    )
    # Create new endpoints for testing
    obj.combine_and_create_new_endpoints({"ANY24": ["DP24", "SM24", "JAW"],
                                          "ANY": ["NC24", "DP24", "SM24", "JAW"]})
    obj.set_well_to_na(endpoint_name="DNC", endpoint_value=1, except_endpoint=["ANY24"])
    obj.remove_endpoints("DNC")
    obj.filter_min_concentration(apply=True)
    return obj

# Test binary fit 
def test_binary_model_fitting(capsys):
    obj = _binary()
    obj.filter_correlation_score(apply = True)
    obj.fit_models(diagnostic_mode = True)
    assert not obj.bmds.empty
    assert obj.bmds_filtered is not None
    assert "fitting models for" in capsys.readouterr().out

# Test where fitting builds plate groups and runs on filtered object
def test_binary_model_fitting_build_plates(capsys):
    obj = _binary()
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    assert not obj.bmds.empty

# Test fit without a correlation filter
def test_binary_model_fitting_no_filtering(capsys):
    obj = _binary()
    obj.filter_correlation_score(apply=False)
    obj.fit_models()
    assert obj.bmds_filtered is not None

# Test fit with invalid parameters
def test_fit_input_validation(capsys):
    obj = _binary()
    obj.filter_correlation_score(apply = False)
    obj.fit_models(gof_threshold = 5, model_selection = "bogus")
    out = capsys.readouterr().out
    assert "gof_threshold must be larger than 0" in out
    assert "only 'lowest BMDL' is supported" in out

# Test proportional fit with removed endpoints using real data
def test_proportional_model_fitting(tmp_path):
        tox = pd.read_csv("data/ToxExample_Long.txt", sep="\t")
        # Keep endpoints with enough concentrations to fit
        tox = tox[tox["endpoint"].isin(["Endpoint1", "202992"])]
        obj = ProportionalClass.ProportionalClass(
            df = tox,
            chemical = "chemical.id",
            concentration = "concentration",
            endpoint = "endpoint",
            response = "response"  
        )
        obj.filter_correlation_score(apply = True)
        obj.fit_models()
        assert obj.bmds is not None

# Test response curve for all models and validation
def test_response_curve():
    obj = _binary()
    obj.filter_correlation_score(apply = False)
    obj.fit_models()
    
    # invalid chemical / endpoint / model
    with pytest.raises(ValueError):
        obj.response_curve(chemical_name = "nope", endpoint_name = "ANY", model = "logistic")
    with pytest.raises(ValueError):
        obj.response_curve(chemical_name = "1", endpoint_name = "nope", model = "logistic")
    with pytest.raises(ValueError):
        obj.response_curve(chemical_name = "1", endpoint_name = "ANY", model = "not_a_model")

    # Every model option
    for model in ["logistic", "gamma", "weibull", "log logistic", "probit", "log probit", "multistage2",  "quantal linear"]:
        obj.response_curve(chemical_name = "1", endpoint_name = "ANY", model = model)

# Test fitting normal and unfit endpoints table
def test_fitting_endpoints_table(tmp_path):
    obj = _binary()
    obj.filter_correlation_score(apply = True)
    obj.fit_models()
    obj.output_fits_table(path=str(tmp_path / "fits.csv"))
    assert not obj.output_res_fits_table.empty

# Test to calculate BMD from an unrecognized model
def test_calculate_bmd_unrecognized_model(capsys):
    result = Calculate_BMD(Model = "not_a_model", params = [1, 2])
    assert result is None
    assert "was not recognized" in capsys.readouterr().out

# Test to calculate BMDL from an unrecognized model
def test_calculate_bmdl_unrecognized_model(capsys):
    obj = _binary()
    obj.filter_correlation_score(apply = False)
    obj.fit_models()
    Data = obj.plate_groups[obj.plate_groups["bmdrc.Endpoint.ID"] == "2 ANY"]
    Data = Data[["concentration", "bmdrc.num.affected", "bmdrc.num.nonna"]].astype(float)

    #using from bmdrc.model_fitting import Logistic
    real_fit = Logistic(Data)

    with pytest.raises(AttributeError):
        Calculate_BMDL("concentration", "not_a_model", real_fit, Data, 1.0, [1, 2])
    assert "was not recognized" in capsys.readouterr().out

# Test to make all models fail using strict gof_threshold
def test_all_models_fail_gof_threshold():
    obj = _binary()
    obj.remove_endpoints("DP24")
    obj.filter_correlation_score(apply=False)
    obj.fit_models(gof_threshold=0.99999)
    assert hasattr(obj, "failed_pvalue_test")

# Test to make all models fail aic_threshold
def test_all_models_fail_aic_threshold():
    obj = _binary()
    obj.filter_correlation_score(apply=False)
    obj.fit_models(gof_threshold=0.0, aic_threshold=0)
    assert hasattr(obj, "failed_pvalue_test")

# Test to calculate bmdl model exception branch
def test_calculate_bmdl_model_exception_branch():
    obj = _binary()
    obj.filter_correlation_score(apply=False)
    obj.fit_models()
    Data=obj.plate_groups[obj.plate_groups["bmdrc.Endpoint.ID"] == "2 ANY"]
    Data=Data[["concentration", "bmdrc.num.affected", "bmdrc.num.nonna"]].astype(float)
    real_fit = Logistic(Data)

    for model in["Logistic", "Gamma", "Weibull", "Log Logistic", "Probit", "Log Probit", 
                 "Multistage2", "Quantal Linear"]:
         result = Calculate_BMDL("concentration", model, real_fit, Data, 1.0, [])
         assert np.isnan(result)