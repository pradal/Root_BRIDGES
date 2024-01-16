from dataclasses import dataclass, asdict

from coupling_demo.model_1.model import Model_1, Model_1_parameters
from coupling_demo.model_2.model import Model_2, Model_2_parameters


# During coupling, we can add new parameters from class inheritance
@dataclass
class Model_1_coupled_parameters(Model_1_parameters):
    new_parameter: float = 3

# Unchanged
@dataclass
class Model_2_coupled_parameters(Model_2_parameters):
    unused: float = 0


class Model_1_Coupled(Model_1):
    """
    This class inherits an existing model to easily edit its initialization, processes and pool balance strategies
    (update or replace)
    """
    def __init__(self, g):
        super().__init__(g)
        # Note : initialized property has to be mutable for proper links between models.
        self.g.update({"new_pool": {1: 1}})
        self.new_pool = self.g["new_pool"]
        self.inputs.update({
            # Common
            "carbon": ["external_property",
                       "external_property_2"]})


    def process_M_1(self, pM1, new_parameter, **kwargs):
        return pM1 + new_parameter

    def process_M_2(self, new_parameter, **kwargs):
        return 10 * new_parameter + self.external_property.link[1]

    def update_pool_2(self):
        self.new_pool[1] += self.process_M_1_value - self.process_M_2_value


# Unchanged
class Model_2_Coupled(Model_2):
    def __init__(self, g):
        super().__init__(g)
