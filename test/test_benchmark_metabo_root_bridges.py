# Public packages
import os, sys, time
import pandas as pd
import multiprocessing as mp
# Model packages
from root_bridges.root_bridges_metabo import Model
# Utility packages
from log.logging import Logger
from analyze.analyze import analyze_data
from initialize import MakeScenarios as ms


def single_run(scenario, outputs_dirpath="outputs", simulation_length=2500, echo=True, log_settings={}):
    root_bridges_metabo = Model(3600, **scenario)

    logger = Logger(model_instance=root_bridges_metabo, outputs_dirpath=outputs_dirpath, 
                    time_step_in_hours=1, logging_period_in_hours=24,
                    echo=echo, **log_settings)
    
    global_performance = pd.DataFrame()
    
    try:
        for _ in range(simulation_length):
            # Placed here also to capture mtg initialization
            t1 = time.time()
            logger()
            t2 = time.time()
            logger.run_and_monitor_model_step()
            t3 = time.time()
            #root_bridges_metabo.run()
            global_performance = pd.concat([global_performance, pd.DataFrame({"time": [root_bridges_metabo.time], "log_time": [t2 - t1], "run_time": [t3 - t2], "total_time": [t3 - t1]})])

    except (ZeroDivisionError, KeyboardInterrupt):
        logger.exceptions.append(sys.exc_info())

    finally:
        logger.stop()
        global_performance.to_excel(os.path.join(outputs_dirpath, "global_performance.xlsx"))
        #analyze_data(scenarios=[os.path.basename(outputs_dirpath)], outputs_dirpath=outputs_dirpath, target_properties=None, **log_settings)


def simulate_scenarios(scenarios, simulation_length=2500, echo=True, log_settings={}):
    processes = []
    max_processes = mp.cpu_count()
    for scenario_name, scenario in scenarios.items():
        
        while len(processes) == max_processes:
            for proc in processes:
                if not proc.is_alive():
                    processes.remove(proc)
            time.sleep(1)

        print(f"[INFO] Launching scenario {scenario_name}...")
        p = mp.Process(target=single_run, kwargs=dict(scenario=scenario, 
                                                      outputs_dirpath=os.path.join("outputs", str(scenario_name)),
                                                      simulation_length=simulation_length,
                                                      echo=echo,
                                                      log_settings=log_settings))
        p.start()
        processes.append(p)


if __name__ == '__main__':
    scenarios = ms.from_table(file_path="inputs/Scenarios_24_06.xlsx", which=["Benchmark_All_Metabo"])
    #simulate_scenarios(scenarios, simulation_length=1, log_settings=Logger.light_log)
    for scenario_name in scenarios:
        scenario = scenarios[scenario_name]
        single_run(
            scenario=scenario, 
            outputs_dirpath=os.path.join("outputs", str(scenario_name)),
            )
