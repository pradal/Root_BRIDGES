from dataclasses import dataclass


# PARAMETERS
@dataclass
class Model_1_parameters:
    pM1: float = 1
    pT1: float = 2


class Model_1:
    def __init__(self, g):
        self.g = g
        self.props = self.g.properties()
        # Note : initialized property has to be mutable for proper links between models.
        self.props.update({"pool_1": {1: 1}})
        self.pool_1 = self.props["pool_1"]
        self.inputs = {}

    # TRANSPORT METHODS
    def process_T_1(self, pT1, **kwargs):
        return pT1

    # METABOLISM METHODS
    def process_M_1(self, pM1, **kwargs):
        return 2*pM1

    # UPDATE POOLS
    def update_pool_1(self):
        self.pool_1[1] += self.process_T_1_value + self.process_M_1_value

    def exchanges_and_balance(self, **kwargs):
        """
        The only expectation for this function to differentiate pools actualization methods with the "update mention"
        """
        process_methods = [func for func in dir(self) if (callable(getattr(self, func)) and '__' not in func and 'process' in func)]
        update_methods = [func for func in dir(self) if (callable(getattr(self, func)) and '__' not in func and 'update' in func)]

        # Compute flows
        for method in process_methods:
            setattr(self, method + "_value", getattr(self, method)(**kwargs))

        # Actualize pools
        for method in update_methods:
            getattr(self, method)()
