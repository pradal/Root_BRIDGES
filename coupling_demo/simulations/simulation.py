from dataclasses import asdict
import pickle
from coupling_demo.coupled_model.coupled_model import Model_1_Coupled, Model_1_coupled_parameters, Model_2_Coupled, Model_2_coupled_parameters
from coupling_demo.coupled_model import converter

# Fake MTG dict
with open("root03360.pckl", "rb") as f:
    g = pickle(f)


# Initialization of the extended models
m1 = Model_1_Coupled(g)
m2 = Model_2_Coupled(g)

# Linking communication variables with different names
converter.link_mtg(receiver=m1, applier=m2, category="carbon", translator=converter.shared_states, same_names=False)

# Computation loop
for k in range(10):
    m1.exchanges_and_balance(**asdict(Model_1_coupled_parameters()))
    m2.exchanges_and_balance(**asdict(Model_2_coupled_parameters()))
    print(m1.external_property.link)
    print(m1.external_property_2.link)
