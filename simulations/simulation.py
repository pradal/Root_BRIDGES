# Public packages
import os, sys, time
import multiprocessing as mp
# Model packages
from root_bridges.root_bridges import Model
# Utility packages
from log.logging import Logger
from analyze.analyze import analyze_data
from initialize.initialize import MakeScenarios as ms


def single_run(scenario, outputs_dirpath="outputs"):
    root_bridges = Model(time_step=3600, **scenario)
    
    logger = Logger(model_instance=root_bridges, outputs_dirpath=outputs_dirpath, 
                    time_step_in_hours=1,
                    logging_period_in_hours=3,
                    recording_images=True, recording_off_screen=True, plotted_property="C_hexose_root", show_soil=True,
                    recording_mtg=False,
                    recording_raw=False,
                    recording_sums=True,
                    recording_performance=True,
                    echo=True)
    
    try:
        for _ in range(2500):
            # Placed here also to capture mtg initialization
            logger()
            logger.run_and_monitor_model_step()
            #root_bridges.run()

    except (ZeroDivisionError, KeyboardInterrupt):
        logger.exceptions.append(sys.exc_info())

    finally:
        logger.stop()
        analyze_data(outputs_dirpath=outputs_dirpath,
                     on_sums=True,
                     on_performance=True,
                     animate_raw_logs=False,
                     target_properties=None
                     )


def simulate_scenarios():
    scenarios = ms.from_table(file_path="inputs/Scenarios_24_05.xlsx", which=["Fischer_1"])
    processes = []
    max_processes = mp.cpu_count()
    for scenario_name, scenario in scenarios.items():
        while len(processes) == max_processes:
            for proc in processes:
                if not proc.is_alive():
                    processes.remove(proc)
            time.sleep(1)

        print(f"[INFO] Launching scenario {scenario_name}...")
        p = mp.Process(target=single_run, kwargs=dict(scenario=scenario, outputs_dirpath=os.path.join("outputs", str(scenario_name))))
        p.start()
        processes.append(p)


if __name__ == '__main__':
    simulate_scenarios()
    # In the end put the system to sleep, windows only
    #os.system("rundll32.exe powrprof.dll,SetSuspendState 0,1,0")
