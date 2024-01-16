from dataclasses import dataclass, asdict

# from coupling_demo.model_1.model import Model_1, Model_1_parameters
from root_cynaps.root_nitrogen import RootNitrogenModel
from root_cynaps.root_water import RootWaterModel

from rhizodep.root_carbon import RootCarbonModel
from rhizodep.root_growth import RootGrowthModel
from rhizodep.root_anatomy import RootAnatomy
from rhizodep.rhizo_soil import SoilModel


# EDITED MODELS


class SoilModelCoupled(SoilModel):
    """
    This class inherits an existing model to easily edit its initialization, processes and pool balance strategies
    (update or replace)
    """
    def __init__(self, g, time_step, **scenario):
        super().__init__(g, time_step, **scenario)


class RootCarbonModelCoupled(RootCarbonModel):
    """
    This class inherits an existing model to easily edit its initialization, processes and pool balance strategies
    (update or replace)
    """
    # --- INPUTS STATE VARIABLES FROM OTHER COMPONENTS : default values are provided if not superimposed by model coupling ---
    # FROM NITROGEN MODEL
    def __init__(self, g, time_step, **scenario):
        super().__init__(g, time_step, **scenario)


class RootNitrogenModelCoupled(RootNitrogenModel):
    """
    This class inherits an existing model to easily edit its initialization, processes and pool balance strategies
    (update or replace)
    """
    def __init__(self, g, time_step, **scenario):
        super().__init__(g, time_step, **scenario)


class RootWaterModelCoupled(RootWaterModel):
    """
    This class inherits an existing model to easily edit its initialization, processes and pool balance strategies
    (update or replace)
    """
    def __init__(self, g, time_step, **scenario):
        super().__init__(g, time_step, **scenario)


# This class inherits an existing model to easily edit its initialization, processes and pool balance strategies (update or replace)
class RootGrowthModelCoupled(RootGrowthModel):
    """
    Doc replacer
    """
    def __init__(self, time_step, **scenario):
        super().__init__(time_step, **scenario)


class RootAnatomyCoupled(RootAnatomy):
    """
    This class inherits an existing model to easily edit its initialization, processes and pool balance strategies
    (update or replace)
    """
    def __init__(self, g, time_step, **scenario):
        super().__init__(g, time_step, **scenario)
