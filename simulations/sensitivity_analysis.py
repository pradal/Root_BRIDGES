# Public packages
import os, sys, time
import numpy as np
import multiprocessing as mp
from SALib.analyze import sobol
# Model packages
from root_bridges.root_bridges import Model
# Utility packages
from log.logging import Logger
from analyze.analyze import analyze_data
from initialize.initialize import MakeScenarios as ms
from initialize.initialize import read_table
import matplotlib.pyplot as plt
from simulations.simulation import simulate_scenarios


def sobol_analysis(problem, results_dirpath, scenarios_names, times=[], outputs=[]):
    result_tables = [read_table(os.join.path(results_dirpath, name, "MTG_properties/MTG_properties_summed/plant_scale_properties.csv")) 
                     for name in scenarios_names]
    
    analyses = {}
    for t in times:
        results = {output: np.array([table[output][t] for table in result_tables]) for output in outputs}
        analyses[t] = {output: sobol.analyze(problem, results[output]) for output in outputs}

    return analyses


if __name__ == '__main__':
    problem, scenarios_filename, scenarios_names = ms.from_factorial_plan("inputs/Factorial_plan_SA.xlsx")
    scenarios = ms.from_table(scenarios_filename, which=scenarios_names)
    simulate_scenarios(scenarios, simulation_length=25, echo=False)
    times=[20]
    outputs=["total_struct_mass"]
    analyses = sobol_analysis(problem=problem, results_dirpath="outputs", scenarios_names=scenarios_names, times=times, outputs=outputs)

    fig, ax = plt.subplots(len(times), len(outputs), figsize=(16, 6))
    for t in range(len(times)):
        for o in range(len(outputs)):
            indices = analyses[times[t]][outputs[o]]
            ax[t][o].bar(problem["names"], indices["S1"], yerr=indices["S1_conf"])
            ax[t][o].set_title(f"t = {times[t]}, {outputs[o]}")
    plt.show()

