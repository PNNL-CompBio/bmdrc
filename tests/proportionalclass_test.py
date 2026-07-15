import pandas as pd
import pytest
import bmdrc
from bmdrc import ProportionalClass

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=SyntaxWarning)

## ProportionalClass Tests ##

# Helper function to create valid proportional data
def make_prop_df():
    return pd.DataFrame({
        "chemical.id": ["Chem1", "Chem1", "Chem1", "Chem2"],
        "concentration": [0, 1, 2, 0],
        "endpoint": ["EP1", "EP1", "EP1", "EP2"],
        "response": [0.0, 0.5, 1.0, 0.25]
    })


# Helper function to create a valid ProportionalClass object
def make_prop_test():
    return ProportionalClass.ProportionalClass(
        df=make_prop_df(),
        chemical="chemical.id",
        concentration="concentration",
        endpoint="endpoint",
        response="response"
    )


# Test that valid proportional data runs without error
def test_proportional_class_valid_init():
    Prop_Test = make_prop_test()

    # Ensure the object was created correctly
    assert isinstance(Prop_Test, ProportionalClass.ProportionalClass)

    # ProportionalClass should add a placeholder plate column
    assert "plate" in Prop_Test.df.columns
    assert set(Prop_Test.df["plate"]) == {"NoPlate"}

    # Ensure inputs were stored correctly
    assert Prop_Test.chemical == "chemical.id"
    assert Prop_Test.concentration == "concentration"
    assert Prop_Test.endpoint == "endpoint"
    assert Prop_Test.response == "response"


# Test wrong inputs for df
def test_proportional_df_errors():
    # df must be a pandas DataFrame
    with pytest.raises(Exception, match="df must be a pandas DataFrame."):
        ProportionalClass.ProportionalClass(
            df="not_a_dataframe",
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )

    # df cannot be empty
    with pytest.raises(Exception, match="df cannot be empty. Please provide a pandas DataFrame."):
        ProportionalClass.ProportionalClass(
            df=pd.DataFrame(),
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )


# Test wrong inputs for chemical
def test_proportional_chemical_errors():
    # chemical must be a string
    with pytest.raises(Exception, match="chemical must be a name of a column in df."):
        ProportionalClass.ProportionalClass(
            df=make_prop_df(),
            chemical=3,
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )

    # chemical must exist in df columns
    with pytest.raises(Exception, match="missing_chemical is not in the column names of df."):
        ProportionalClass.ProportionalClass(
            df=make_prop_df(),
            chemical="missing_chemical",
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )

    # chemical cannot be an unacceptable name
    with pytest.raises(Exception, match="bmdrc.Well.ID is not a permitted name. Please rename this column."):
        df2 = make_prop_df().rename({"chemical.id": "bmdrc.Well.ID"}, axis=1)
        ProportionalClass.ProportionalClass(
            df=df2,
            chemical="bmdrc.Well.ID",
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )

    # chemical names cannot contain spaces
    with pytest.raises(Exception, match="spaces are not permitted in chemical names. Check this column in your dataframe."):
        df2 = make_prop_df()
        df2.loc[0, "chemical.id"] = "Chem 1"
        ProportionalClass.ProportionalClass(
            df=df2,
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )


# Test wrong inputs for concentration
def test_proportional_concentration_errors():
    # concentration must be a string
    with pytest.raises(Exception, match="concentration must be a name of a column in df."):
        ProportionalClass.ProportionalClass(
            df=make_prop_df(),
            chemical="chemical.id",
            concentration=22,
            endpoint="endpoint",
            response="response"
        )

    # concentration must exist in df columns
    with pytest.raises(Exception, match="missing_conc is not in the column names of df."):
        ProportionalClass.ProportionalClass(
            df=make_prop_df(),
            chemical="chemical.id",
            concentration="missing_conc",
            endpoint="endpoint",
            response="response"
        )

    # concentration cannot be an unacceptable name
    with pytest.raises(Exception, match="bmdrc.frac.affected is not a permitted name. Please rename this column."):
        df2 = make_prop_df().rename({"concentration": "bmdrc.frac.affected"}, axis=1)
        ProportionalClass.ProportionalClass(
            df=df2,
            chemical="chemical.id",
            concentration="bmdrc.frac.affected",
            endpoint="endpoint",
            response="response"
        )


# Test wrong inputs for endpoint
def test_proportional_endpoint_errors():
    # endpoint must be a string
    with pytest.raises(Exception, match="endpoint must be a name of a column in df."):
        ProportionalClass.ProportionalClass(
            df=make_prop_df(),
            chemical="chemical.id",
            concentration="concentration",
            endpoint=55,
            response="response"
        )

    # endpoint must exist in df columns
    with pytest.raises(Exception, match="missing_endpoint is not in the column names of df."):
        ProportionalClass.ProportionalClass(
            df=make_prop_df(),
            chemical="chemical.id",
            concentration="concentration",
            endpoint="missing_endpoint",
            response="response"
        )

    # endpoint cannot be an unacceptable name
    with pytest.raises(Exception, match="bmdrc.filter is not a permitted name. Please rename this column."):
        df2 = make_prop_df().rename({"endpoint": "bmdrc.filter"}, axis=1)
        ProportionalClass.ProportionalClass(
            df=df2,
            chemical="chemical.id",
            concentration="concentration",
            endpoint="bmdrc.filter",
            response="response"
        )

    # endpoint names cannot contain spaces
    with pytest.raises(Exception, match="spaces are not permitted in endpoint names. Check this column in your dataframe."):
        df2 = make_prop_df()
        df2.loc[0, "endpoint"] = "EP 1"
        ProportionalClass.ProportionalClass(
            df=df2,
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )


# Test wrong inputs for response
def test_proportional_response_errors():
    # response must be a string
    with pytest.raises(Exception, match="response must be a name of a column in df."):
        ProportionalClass.ProportionalClass(
            df=make_prop_df(),
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response=22
        )

    # response must exist in df columns
    with pytest.raises(Exception, match="missing_response is not in the column name of df."):
        ProportionalClass.ProportionalClass(
            df=make_prop_df(),
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response="missing_response"
        )

    # response cannot be an unacceptable name
    with pytest.raises(Exception, match="bmdrc.filter.reason is not a permitted name. Please rename this column."):
        df2 = make_prop_df().rename({"response": "bmdrc.filter.reason"}, axis=1)
        ProportionalClass.ProportionalClass(
            df=df2,
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response="bmdrc.filter.reason"
        )

    # response values must be between 0 and 1
    with pytest.raises(Exception, match="The response column must range in values from 0 to 1."):
        df2 = make_prop_df()
        df2.loc[0, "response"] = 1.5
        ProportionalClass.ProportionalClass(
            df=df2,
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )

    # response values must not be below 0
    with pytest.raises(Exception, match="The response column must range in values from 0 to 1."):
        df2 = make_prop_df()
        df2.loc[0, "response"] = -0.1
        ProportionalClass.ProportionalClass(
            df=df2,
            chemical="chemical.id",
            concentration="concentration",
            endpoint="endpoint",
            response="response"
        )


# Test wrapper method: remove_endpoints
def test_proportional_wrapper_methods_call_imported_functions(monkeypatch):
    Prop_Test = make_prop_test()
    called = []

    def fake_remove_endpoints(self, endpoint_name):
        called.append(("remove_endpoints", endpoint_name))

    def fake_min_concentration(self, count, apply, diagnostic_plot):
        called.append(("min_concentration", count, apply, diagnostic_plot))

    def fake_correlation_score(self, score, apply, diagnostic_plot, direction):
        called.append(("correlation_score", score, apply, diagnostic_plot, direction))

    def fake_fit_the_models(self, gof_threshold, aic_threshold, model_selection, diagnostic_mode):
        called.append(("fit_the_models", gof_threshold, aic_threshold, model_selection, diagnostic_mode))

    def fake_gen_response_curve(self, chemical_name, endpoint_name, model, steps):
        called.append(("gen_response_curve", chemical_name, endpoint_name, model, steps))

    def fake_benchmark_dose(self, path):
        called.append(("benchmark_dose", path))

    def fake_dose_table(self, path):
        called.append(("dose_table", path))

    def fake_fits_table(self, path):
        called.append(("fits_table", path))

    def fake_report_binary(self, out_folder, report_name, file_type):
        called.append(("report_binary", out_folder, report_name, file_type))

    monkeypatch.setattr(ProportionalClass, "remove_endpoints", fake_remove_endpoints)
    monkeypatch.setattr(ProportionalClass, "min_concentration", fake_min_concentration)
    monkeypatch.setattr(ProportionalClass, "correlation_score", fake_correlation_score)
    monkeypatch.setattr(ProportionalClass, "fit_the_models", fake_fit_the_models)
    monkeypatch.setattr(ProportionalClass, "gen_response_curve", fake_gen_response_curve)
    monkeypatch.setattr(ProportionalClass, "benchmark_dose", fake_benchmark_dose)
    monkeypatch.setattr(ProportionalClass, "dose_table", fake_dose_table)
    monkeypatch.setattr(ProportionalClass, "fits_table", fake_fits_table)
    monkeypatch.setattr(ProportionalClass, "report_binary", fake_report_binary)

    # Ensure the preprocessing wrapper delegates to the imported helper.
    Prop_Test.remove_endpoints(["EP1"])

    # Ensure the filtering wrappers delegate with the provided arguments.
    Prop_Test.filter_min_concentration(count = 2, apply = True, diagnostic_plot = True)
    Prop_Test.filter_correlation_score(score = 0.5, apply = True, diagnostic_plot = False, direction = "above")

    # Ensure the model-fitting wrappers delegate to the imported fitting helpers.
    Prop_Test.fit_models(gof_threshold = 0.2, aic_threshold = 3, model_selection = "lowest BMDL", diagnostic_mode = True)
    Prop_Test.response_curve("Chem1", "EP1", "logistic", 12)

    # Ensure the output wrappers delegate to the imported output helpers.
    Prop_Test.output_benchmark_dose("benchmark.csv")
    Prop_Test.output_dose_table("dose.csv")
    Prop_Test.output_fits_table("fits.csv")
    Prop_Test.report("reports", "Dose Report", ".html")

    assert called == [
        ("remove_endpoints", ["EP1"]),
        ("min_concentration", 2, True, True),
        ("correlation_score", 0.5, True, False, "above"),
        ("fit_the_models", 0.2, 3, "lowest BMDL", True),
        ("gen_response_curve", "Chem1", "EP1", "logistic", 12),
        ("benchmark_dose", "benchmark.csv"),
        ("dose_table", "dose.csv"),
        ("fits_table", "fits.csv"),
        ("report_binary", "reports", "Dose Report", ".html"),
    ]