from root_bridges.root_growth import RootGrowthModelCoupled

# Utilities
from metafspm.composite_wrapper import CompositeModel
from metafspm.component_factory import Choregrapher


class Model(CompositeModel):
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

        # DECLARE GLOBAL SIMULATION TIME STEP
        Choregrapher().add_simulation_time_step(time_step)
        self.time = 0
        parameters = scenario["parameters"]["root_bridges"]
        self.input_tables = scenario["input_tables"]

        # INIT INDIVIDUAL MODULES
        if len(scenario["input_mtg"]) > 0:
            self.root_growth = RootGrowthModelCoupled(scenario["input_mtg"]["root_mtg_file"], time_step, **parameters)
        else:
            self.root_growth = RootGrowthModelCoupled(time_step, **parameters)
        self.g = self.root_growth.g

        # EXPECTED !
        self.models = (self.root_growth,)
        self.data_structures = {"root": self.g}


    def run(self):
        self.apply_input_tables(tables=self.input_tables, to=self.models, when=self.time)

        # Compute root growth from resulting states
        self.root_growth()

        self.time += 1

