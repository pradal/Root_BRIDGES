from dataclasses import dataclass


# PARAMETERS
@dataclass
class Model_2_parameters:
    p1: float = 1


class Model_2:
    def __init__(self, g):
        self.g = g
        # Note : initialized property has to be mutable for proper links between models.
        self.g.update({"pool_3": {1: 1}, "pool_4": {1: 1}})
        self.pool_3 = self.g["pool_3"]
        self.pool_4 = self.g["pool_4"]
        self.inputs = {}

    def process_1(self, p1, **kwargs):
        return 5*p1

    # UPDATE POOLS
    def update_pool_1(self):
        self.pool_3[1] += self.process_1_value
        self.pool_4[1] -= self.process_1_value ** 0.5

    def exchanges_and_balance(self, **kwargs):
        """
        The only expectation for this function to differentiate pools actualization methods with the "update mention"
        """
        process_methods = [func for func in dir(self) if (callable(getattr(self, func)) and '__' not in func and 'update' not in func and func != 'exchanges_and_balance')]
        update_methods = [func for func in dir(self) if (callable(getattr(self, func)) and '__' not in func and 'update' in func)]

        # Compute flows
        for method in process_methods:
            setattr(self, method + "_value", getattr(self, method)(**kwargs))

        # Actualize pools
        for method in update_methods:
            getattr(self, method)()
