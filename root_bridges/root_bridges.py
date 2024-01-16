import os
import pickle

from root_bridges.coupled_models import (
    # MODELS
    RootNitrogenModelCoupled,
    RootCarbonModelCoupled,
    RootWaterModelCoupled,
    RootGrowthModelCoupled,
    RootAnatomyCoupled,
)

from rhizodep.rhizo_soil import SoilModel
from Data_enforcer.shoot import ShootModel

from root_cynaps.wrapper import ModelWrapper


class Model(ModelWrapper):
    """
    Root-BRIDGES model

    Use guideline :
    1. store in a variable Model(g, time_step) to initialize the model, g being an openalea.MTG() object and time_step an time interval in seconds.

    2. print Model.documentation for more information about editable model parameters (optional).

    3. Use Model.scenario(**dict) to pass a set of scenario-specific parameters to the model (optional).

    4. Use Model.run() in a for loop to perform the computations of a time step on the passed MTG File
    """

    def __init__(self, time_step: int, **scenario):
        """
        DESCRIPTION
        ----------
        __init__ method of the model. Initializes the thematic modules and link them.

        :param g: the openalea.MTG() instance that will be worked on. It must be representative of a root architecture.
        :param time_step: the resolution time_step of the model in seconds.
        """

        # INIT INDIVIDUAL MODULES
        self.root_growth = RootGrowthModelCoupled(time_step, **scenario)
        self.g = self.root_growth.g
        self.root_anatomy = RootAnatomyCoupled(self.g, time_step, **scenario)
        self.root_carbon = RootCarbonModelCoupled(self.g, time_step, **scenario)
        self.root_water = RootWaterModelCoupled(self.g, time_step)
        self.root_nitrogen = RootNitrogenModelCoupled(self.g, time_step)
        self.soil = SoilModel(self.g, time_step, **scenario)

        # Initialisation of Shoot modules
        self.shoot = ShootModel(self.g)

        self.models = (self.root_growth, self.root_anatomy, self.root_carbon, self.root_water, self.root_nitrogen, self.soil)

        # LINKING MODULES
        # Get or build translator matrix
        if not os.path.isfile("translator.pckl"):
            print("NOTE : You will now have to provide information about shared variables between the modules composing this model :\n")
            self.translator_matrix_builder()
        with open("translator.pckl", "rb") as f:
            translator = pickle.load(f)

        # Actually link modules together
        self.link_around_mtg(translator)

        # Some initialization must be performed AFTER linking modules
        self.soil.post_coupling_init()
        self.root_water.post_coupling_init()
        self.root_carbon.post_coupling_init()
        self.root_nitrogen.post_coupling_init()
        self.root_growth.post_coupling_init()
        self.root_anatomy.post_coupling_init()

        # Step number to read datatable for data enforcing
        self.step = 1

    def run(self):
        # Update environment boundary conditions
        self.soil.run_exchanges_and_balance()

        # Compute shoot flows and state balance
        self.shoot.run_exchanges_and_balance(time=self.step)

        # Compute state variations for water and then carbon and nitrogen
        self.root_water.run_exchanges_and_balance()
        self.root_carbon.run_exchanges_and_balance()
        self.root_nitrogen.run_exchanges_and_balance()

        # Compute root growth from resulting states
        self.root_growth.run_time_step_growth()
        # Update topological surfaces and volumes based on other evolved structural properties
        self.root_anatomy.run_actualize_anatomy()

        self.step += 1
