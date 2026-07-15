import numpy as np
import pandas as pd
import pytest
import bmdrc
from bmdrc import ContinuousClass
from contextlib import contextmanager

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)


# Test dataset for valid input data
def test_dataset_ContinuousClass():
    return pd.read_csv("data/Continuous.txt", sep="\t").drop("Notes", axis=1)

def _make(**overrides):
    """ Build a ContinuousClass instance with default parameters, allowing overrides. """
    args = dict(df=test_dataset_ContinuousClass(), chemical="Chemical ID", concentration="Concentration_uM", endpoint="Endpoint", response="Measurement")
    args.update(overrides)
    return ContinuousClass.ContinuousClass(**args)

# ContinuousClass Tests construction 
def test_cotinuous_construct():
    obj = _make()
    assert isinstance(obj, ContinuousClass.ContinuousClass)

## Test dataframe is a pandas DataFrame
def test_must_be_dataframe():
    with pytest.raises(Exception, match="df must be a pandas DataFrame"):
        _make(df="not a dataframe")

# Test dataframe is not empty
def test_dataframe_not_empty():
    with pytest.raises(Exception, match="df cannot be empty. Please provide a pandas DataFrame."):
        _make(df=pd.DataFrame())

# Test chemical column is a string
def test_chemical_column_is_string():
    with pytest.raises(Exception, match="chemical must be a name of a column in df."):
        _make(chemical=3)

# Test chemical column exists in dataframe
def test_chemical_column_exists():
    with pytest.raises(Exception, match="bmrdc.filter is not in the column names of df."):
        _make(chemical="bmrdc.filter")

# Test chemical.id rejects reserved name
def test_chemical_column_name():
    df = test_dataset_ContinuousClass().rename({"Chemical ID": "bmrdc.filter"}, axis=1)
    with pytest.raises(Exception, match="bmrdc.filter is not in the column names of df."):
        _make(df = df, chemical=" bmrdc.filter")

def test_chemical_column_rejects_reserved_name():
    df = test_dataset_ContinuousClass().rename({"Chemical ID": "bmdrc.Well.ID"}, axis=1)
    with pytest.raises(Exception, match="bmdrc.Well.ID is not a permitted name. Please rename this column."):
        _make(df=df, chemical="bmdrc.Well.ID")

# Test chemical rejects spaces
def test_chemical_column_no_spaces():
    df = test_dataset_ContinuousClass().copy()
    df["Chemical ID"] = df["Chemical ID"].astype(str)
    df.loc[df.index[0], "Chemical ID"] = " 1"
    with pytest.raises(Exception, match="spaces are not permitted in chemical names. Check this column in your dataframe."):
        _make(df = df)

# Test concentration column is a string
def test_concentration_column_is_string():
    with pytest.raises(Exception, match="concentration must be a name of a column in df."):
        _make(concentration=3)

# Test concentration column exists in dataframe
def test_concentration_column_exists():
    with pytest.raises(Exception, match="bmrdc.filter is not in the column names of df."):
        _make(concentration="bmrdc.filter")

def test_concentration_column_rejects_reserved_name():
    df = test_dataset_ContinuousClass().rename({"Concentration_uM": "bmdrc.num.nonna"}, axis=1)
    with pytest.raises(Exception, match="bmdrc.num.nonna is not a permitted name. Please rename this column."):
        _make(df=df, concentration="bmdrc.num.nonna")

# Test endpoint value is a string
def test_endpoint_column_is_string():
    with pytest.raises(Exception, match="endpoint must be a name of a column in df."):
        _make(endpoint=77)

# Test endpoint column exists in dataframe
def test_endpoint_column_exists():
    with pytest.raises(Exception, match="bmrdc.filter is not in the column names of df."):
        _make(endpoint="bmrdc.filter")

def test_endpoint_column_rejects_reserved_name():
    df = test_dataset_ContinuousClass().rename({"Endpoint": "bmdrc.num.tot"}, axis=1)
    with pytest.raises(Exception, match="bmdrc.num.tot is not a permitted name. Please rename this column."):
        _make(df=df, endpoint="bmdrc.num.tot")

# Test enpoint column rejects spaces
def test_endpoint_column_no_spaces():
    df = test_dataset_ContinuousClass()
    df.loc[df.index[0], "Endpoint"] = " Endpoint6"
    with pytest.raises(Exception, match="spaces are not permitted in endpoint names. Check this column in your dataframe."):
        _make(df = df)

# Test response must be a string
def test_response_column_is_string():
    with pytest.raises(Exception, match="response must be a name of a column in df."):
        _make(response=33)

# Test response column exists in dataframe
def test_response_column_exists():
    with pytest.raises(Exception, match="bmrdc.filter is not in the column name of df."):
        _make(response="bmrdc.filter")

# Test response column rejects spaces
def test_response_column_no_spaces():
    df = test_dataset_ContinuousClass().rename({"Measurement": " measurement value"}, axis=1)
    with pytest.raises(Exception, match="bmrdc.filter is not in the column name of df."):
        _make(df = df, response="bmrdc.filter")

def test_response_column_rejects_reserved_name():
    df = test_dataset_ContinuousClass().rename({"Measurement": "bmdrc.num.affected"}, axis=1)
    with pytest.raises(Exception, match="bmdrc.num.affected is not a permitted name. Please rename this column."):
        _make(df=df, response="bmdrc.num.affected")

@contextmanager
def stub_helper(name, stub_function):
    """ Context manager to temporarily replace a method in ContinuousClass with a stub function. """
    original_method = getattr(ContinuousClass, name)
    setattr(ContinuousClass, name, stub_function)
    try:
        yield
    finally:
        setattr(ContinuousClass, name, original_method)

# Test remove endpoints using helper
def test_remove_endpoints():
    calls = []
    with stub_helper("remove_endpoints", 
                     lambda obj, endpoints_name: calls.append(endpoints_name)):
        _make().remove_endpoints("Endpoint1")
    assert calls == ["Endpoint1"]

# Test filter min concentration using helper
def test_filter_min_concentration():
    calls = []
    with stub_helper("min_concentration", 
                     lambda obj, count, apply, diagnostic_plot:
                     calls.append((count,apply, diagnostic_plot))):
        _make().filter_min_concentration(count=3, apply=True, diagnostic_plot=False)
    assert calls == [(3, True, False)]

# Test filter correlation score using helper
def test_filter_correlation_score():
    calls = []
    with stub_helper("correlation_score", 
                     lambda obj, score, apply, diagnostic_plot, direction:
                     calls.append((score, apply, diagnostic_plot, direction))):
        _make().filter_correlation_score(score = 0.2, apply = True, diagnostic_plot = False, direction = "below")
    assert calls == [(0.2, True, False, "below")]

# Test filter negative control using helper
def test_filter_negative_control():
    calls = []
    with stub_helper("negative_control_continuous",
                     lambda obj, apply, diagnostic_plot:
                     calls.append((apply, diagnostic_plot))):
        _make().filter_negative_control(apply= True, diagnostic_plot= False)
    assert calls == [(True, False)]

# Test fit models using helper
def test_fit_models():
    calls = []    
    with stub_helper("fit_continuous_models",
                     lambda obj, fixed_intercept, aic_threshold, model_selection, diagnostic_mode:
                     calls.append((fixed_intercept,aic_threshold, model_selection, diagnostic_mode))):
        _make().fit_models(fixed_intercept=0, aic_threshold=2, model_selection= "lowest BMDL", diagnostic_mode= False)
    assert calls == [(0,2,"lowest BMDL",False)]

# Test response curve using helper
def test_response_curve():
    calls = []
    with stub_helper("gen_response_curve",
                     lambda obj, chemical_name, endpoint_name, model, fixed_intercept, add_bmds, steps:
                     calls.append((chemical_name, endpoint_name, model, fixed_intercept, add_bmds,steps))):
        _make().response_curve(chemical_name="12", endpoint_name="Endpoint1", model="linear", fixed_intercept=0)
    assert calls == [("12", "Endpoint1", "linear", 0, False, 10)]

# Test output benchmark using helper
def test_output_benchmark():
    calls = []
    with stub_helper("benchmark_dose",
                     lambda obj, path:
                     calls.append((path))):
        _make().output_benchmark_dose(path="x.csv")
    assert calls == ["x.csv"]

# Test output dose using helpers 
def test_output_dose():
    calls = []
    with stub_helper("dose_table",
                     lambda obj, path:
                     calls.append((path))):
        _make().output_dose_table(path = "y.csv")
    assert calls == ["y.csv"]

# Test output fits table using helper
def test_output_fits_table():
    calls = []
    with stub_helper("fits_table",
                     lambda obj, fixed_intercept, path:
                     calls.append((fixed_intercept, path))):
        _make().output_fits_table(fixed_intercept=0, path="z.csv")
    assert calls == [(0, "z.csv")]

# Test reports using helpers    
def test_reports():
    calls = []
    with stub_helper("report_binary",
                     lambda obj, out_folder, report_name, file_type:
                     calls.append((out_folder, report_name, file_type))):
        _make().report(out_folder="outdir")
    assert calls == [("outdir", "Benchmark Dose Curves", ".md")]