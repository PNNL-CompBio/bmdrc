import pandas as pd

__author__ = "David Degnan"

def make_groups(self):
    '''
    Assign groups based on the chemical, concentration, plate, and endpoint.
    This step should be completed after all pre-processing steps are finished. 
    '''

    # If there is a bmdrc well id, drop it
    if "bmdrc.Well.ID" in self.df.columns.tolist():
        self.plate_groups = self.df.drop([self.well, "bmdrc.Well.ID"], 1).groupby(by = [self.chemical, self.concentration, self.plate, self.endpoint], as_index = False)
    else:
        self.plate_groups = self.df.drop([self.well], 1).groupby(by = [self.chemical, self.concentration, self.plate, self.endpoint], as_index = False)
        
    # Get the number of samples per group
    num_tot_samples = self.plate_groups.size().rename(columns = {"size": "bmdrc.num.tot"})

    # Get the number of non-na samples per groups
    num_nonna = self.plate_groups.count().rename(columns = {"value": "bmdrc.num.nonna"})

    # Get the number affected
    num_affected = self.plate_groups.sum().rename(columns = {"value": "bmdrc.num.affected"})

    # Merge to create missingness dataframe
    self.plate_groups = pd.merge(pd.merge(num_tot_samples, num_nonna), num_affected)

    # Create IDs of chemical.id, plate.id, and endpoint in plate_groups 
    self.plate_groups["bmdrc.Plate.ID"] = self.plate_groups[self.chemical].astype(str) + " " + \
                                            self.plate_groups[self.plate].astype(str) + " " + \
                                            self.plate_groups[self.endpoint].astype(str)