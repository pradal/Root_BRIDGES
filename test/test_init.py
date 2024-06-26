from root_bridges.root_bridges import Model
from initialize.initialize import MakeScenarios as ms

def test_init():
    scenarios = ms.from_table(file_path="inputs/Scenarios_24_06.xlsx", which=["Reference_Fischer"])
    for name, scenario in scenarios.items():
        model = Model(time_step=3600, scenario=scenario)

if __name__ == "__main__":
    test_init()
