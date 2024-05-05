import os
from root_bridges.root_bridges import Model
from log.logging import Logger
from analyze.analyze import analyze_data
from initialize.initialize import MakeScenarios as ms


def single_run(scenario, outputs_dirpath="outputs"):
    root_bridges = Model(time_step=3600, **scenario)
    
    logger = Logger(model_instance=root_bridges, outputs_dirpath=outputs_dirpath, 
                    time_step_in_hours=1,
                    logging_period_in_hours=20,
                    recording_images=True, plotted_property="C_hexose_root", show_soil=True,
                    recording_mtg=True,
                    recording_raw=True,
                    recording_sums=True,
                    recording_performance=True,
                    echo=True)
    
    for step in range(61):
        # Placed here also to capture mtg initialization
        logger()
        root_bridges.run()
    
    logger.stop()
    analyze_data(outputs_dirpath=outputs_dirpath, 
                 on_sums=True,
                 on_performance=True,
                 animate_raw_logs=True,
                 target_properties=[]
                 )
    
def test_apply_scenarios():
    scenarios = ms.from_table(file_path="inputs/Scenarios.xlsx", which=["T2"])
    for scenario_name, scenario in scenarios.items():
        single_run(scenario=scenario, outputs_dirpath=os.path.join("outputs", str(scenario_name)))


test_apply_scenarios()
