import pandas as pd
import numpy as np 
from bmdrc import BinaryClass
from bmdrc import LPRClass
from bmdrc import ProportionalClass
from bmdrc import ContinuousClass

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)


# binary class model
def _binary(apply_filter=True):
    obj = BinaryClass.BinaryClass(
        df= pd.read_csv("data/Binary_Simplified_Long.csv"),
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
    if apply_filter:
        obj.filter_correlation_score(apply=True)
    obj.fit_models()
    return obj

#Test benchmark with filtered endpoints
def test_benchmark_filtered(tmp_path):
    obj = _binary(apply_filter=True)
    assert obj.bmds_filtered is not None
    obj.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    written = pd.read_csv(str(tmp_path / "bmd.csv"))
    assert not written.empty

#Test benchmark with no filteing
def test_benchmark_no_filter(tmp_path):
    # New obj for no filter
    obj = BinaryClass.BinaryClass(
        df= pd.read_csv("data/Binary_Simplified_Long.csv"),
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
    obj.remove_endpoints(["DNC", "NC24"])
    obj.fit_models()
    assert obj.bmds_filtered is None
    obj.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    assert not pd.read_csv(str(tmp_path / "bmd.csv")).empty

# Test dose table output
def test_dose_table(tmp_path):
    obj = _binary(apply_filter=True)
    obj.output_dose_table(path=str(tmp_path / "dose.csv"))
    assert not pd.read_csv(str(tmp_path / "dose.csv")).empty

# Test markdown report for a binary class object
def test_report_binary_markdown(tmp_path):
    obj = _binary(apply_filter=True)
    obj.report(out_folder=str(tmp_path), report_name="binary_report", file_type=".md")

# Test to ensure report creates output folder if it does not exist
def test_report_creates_missing_folder(tmp_path): 
    obj = _binary(apply_filter=True)
    new_folder = str(tmp_path / "does_not_exist_yet")
    obj.report(out_folder=new_folder, report_name="rpt", file_type=".md")

# Save and process LPR data example
def _lpr():
    obj = LPRClass.LPRClass(
        df = pd.read_csv("data/LPR_Long.csv"),
        chemical="chemical.id",
        plate="plate.id",
        well="well",
        concentration="conc",
        time="variable",
        value="value",        
    )
    obj.fit_models()
    return obj 

# Save and process proportional class
def _proportional():
    df = pd.read_csv("data/ToxExample_Long.txt", sep="\t")
    # keep endpoints with enough concentrations to fit
    df = df[df["endpoint"].isin(["Endpoint1", "202992"])]
    obj = ProportionalClass.ProportionalClass(
        df=df, 
        chemical="chemical.id",
        concentration="concentration",
        endpoint="endpoint",
        response="response"        
    )        
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    return obj

# Test LPR class of the markdown report
def test_report_lpr_markdown(tmp_path):
    _lpr().report(out_folder=str(tmp_path), report_name="lpr_report", file_type=".md")

# Test proportional class of the markdown report
def test_report_proportional_markdwon(tmp_path):
    _proportional().report(out_folder=str(tmp_path), report_name="prop_report", file_type=".md")

# Test LPR benchmark dose and dose table outputs
def test_lpr_outputs(tmp_path):
    obj = _lpr()
    obj.output_benchmark_dose(path=str(tmp_path / "lpr_bmd.csv"))
    obj.output_dose_table(path=str(tmp_path / "lpr_dose.csv"))
    assert not pd.read_csv(str(tmp_path / "lpr_bmd.csv")).empty

# Test json report file type ( .json export branch)
def test_report_json(tmp_path):
    obj = _binary(apply_filter=True)
    obj.report(out_folder=str(tmp_path), report_name="binary.json", file_type=".json")\

# Test report with negative control filter
def test_report_negative_control_filter(tmp_path):
    df = pd.read_csv("data/Binary_Simplified_Long.csv")
    obj = BinaryClass.BinaryClass(
        df=df, 
        chemical="chemical.id",
        plate="plate.id",
        well="well",
        concentration="concentration",
        endpoint="endpoint",
        value="value",
        format="long"
    )
    obj.combine_and_create_new_endpoints({"ANY24": ["DP24", "SM24", "JAW"],
                                          "ANY": ["NC24", "DP24", "SM24", "JAW"]})
    obj.set_well_to_na(endpoint_name="DNC", endpoint_value=1, except_endpoint=["ANY24"])
    obj.remove_endpoints("DNC")
    obj.filter_min_concentration(apply=True)
    obj.filter_negative_control(apply=True, diagnostic_plot=True)
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    obj.report(out_folder=str(tmp_path), report_name="nc_report", file_type=".md")    

# Test benchmark dose rebuilds plate groups
def test_benchmark_dose_rebuilds_plate_groups(tmp_path):
    obj=_binary()
    obj.filter_correlation_score(apply=False)
    obj.fit_models()
    del obj.plate_groups
    obj.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    assert hasattr(obj, "plate_groups")


def test_benchmark_dose_collapses_duplicate_fail_rows(tmp_path):
    obj = _binary(apply_filter=True)
    assert obj.bmds_filtered is not None

    # Force overlap: same endpoint appears in both filtered and GOF-fail tables.
    duplicated_endpoint = obj.bmds_filtered["bmdrc.Endpoint.ID"].iloc[0]
    obj.failed_pvalue_test = [duplicated_endpoint]

    obj.output_benchmark_dose(path=str(tmp_path / "bmd.csv"))
    written = pd.read_csv(str(tmp_path / "bmd.csv"))
    matching_rows = written[written["bmdrc.Endpoint.ID"] == duplicated_endpoint]

    assert len(matching_rows) == 1
    modeled_flag = matching_rows["Modeled_Flag"].iloc[0]
    assert "Fail - GOF check" in modeled_flag
    assert ("Fail - other filter" in modeled_flag) or ("Fail - correlation score filter" in modeled_flag)

# Test to cover functions "this step was not conducted"
def test_report_step_was_not_conducted(tmp_path):
        obj = BinaryClass.BinaryClass(
        df= pd.read_csv("data/Binary_Simplified_Long.csv"),
        chemical="chemical.id",
        plate="plate.id",
        well="well",
        concentration="concentration",
        endpoint="endpoint",
        value="value",
        format="long"
        )
        obj.combine_and_create_new_endpoints({"ANY24": ["DP24", "SM24", "JAW"],
                                          "ANY": ["NC24", "DP24", "SM24", "JAW"]})
        obj.report(out_folder=str(tmp_path), report_name="step_not_conducted", file_type=".md")

# Test to cover full report path (negative control + correlation filter + fit)
def test_report_full_pipeline(tmp_path):
    obj = _binary()
    obj.filter_negative_control(apply=True)
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    obj.report(out_folder=str(tmp_path), report_name="full", file_type=".md")

# Test binary report
def test_report_binary(tmp_path):
    obj = ContinuousClass.ContinuousClass(
        df=pd.read_csv("data/Continuous.txt", sep="\t"),
        chemical="Chemical ID",
        concentration="Concentration_uM",
        endpoint="Endpoint",
        response="Measurement"
    )
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    obj.report(out_folder=str(tmp_path), report_name="cont", file_type=".md")

# Test to cover "else: endpoints + None" branch
def test_report_well_na_no_except(tmp_path):
        obj = BinaryClass.BinaryClass(
        df= pd.read_csv("data/Binary_Simplified_Long.csv"),
        chemical="chemical.id",
        plate="plate.id",
        well="well",
        concentration="concentration",
        endpoint="endpoint",
        value="value",
        format="long"
        )
        obj.combine_and_create_new_endpoints({"ANY24": ["DP24", "SM24", "JAW"],
                                          "ANY": ["NC24", "DP24", "SM24", "JAW"]})
        obj.set_well_to_na(endpoint_name="DNC", endpoint_value=1, except_endpoint=None)
        obj.report(out_folder=str(tmp_path), report_name="r", file_type=".md")

# Covers the negative-control-filter applied branch in report
def test_report_with_negative_control_applied(tmp_path):
    obj = _binary()
    obj.filter_negative_control(apply=True, diagnostic_plot=True)
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    obj.report(out_folder=str(tmp_path), report_name="r", file_type=".md")

# Covers the negative-control-filter removed
def test_report_with_negative_control_removed(tmp_path):
    obj = _binary()
    obj.filter_negative_control(percentage=30, apply=True, diagnostic_plot=True)
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    obj.report(out_folder=str(tmp_path), report_name="r", file_type=".md")

# Covers continous (no-"value") else branch 
def test_report_continuous_no_value_branch(tmp_path):
    obj = ContinuousClass.ContinuousClass(
        df=pd.read_csv("data/Continuous.txt", sep="\t"),
        chemical="Chemical ID",
        concentration="Concentration_uM",
        endpoint="Endpoint",
        response="Measurement"
    )
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    obj.report(out_folder=str(tmp_path), report_name="r", file_type=".md")

# Test report benchmark dose
def test_report_benchmark_dose(tmp_path):
    obj = _binary()
    obj.filter_correlation_score(apply=True)
    obj.fit_models()
    obj.report(out_folder=str(tmp_path), report_name="r", file_type=".md")

