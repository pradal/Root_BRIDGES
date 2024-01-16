import os
import xarray as xr
from time import time

from root_bridges import Model

from tools.mtg_dict_to_xarray import mtg_to_dataset
from root_bridges.output_properties import state_extracts, flow_extracts, global_state_extracts, global_flow_extracts

from statistical_tools.main import launch_analysis


def run_coupled_simulation(z_soil_Nm_max, output_path, time_step=3600, steps_number=500, logging=True, running=True,
                           max_time_steps_for_memory=100):
    # TODO identify missing variables to share between projects
    Loc = locals()
    real_parameters = ["output_path", "time_step", "steps_number", "max_time_steps_for_memory"]
    scenario = dict([(key, value) for key, value in Loc.items() if key not in real_parameters])

    if running:
        # We record the starting time of the simulation:
        t_start = time()

        # Initialization of the model
        root_bridges = Model(time_step=time_step)

        g = root_bridges.g

        if logging:
            # Logging
            log_outputs = {}
            for d in [state_extracts, flow_extracts, global_state_extracts, global_flow_extracts]:
                log_outputs.update(d)
            os.mkdir(output_path)
            time_xrs = [mtg_to_dataset(g, variables=log_outputs, time=0)]


        # Scheduler : actual computation loop
        for step in range(steps_number):
            root_bridges.run()

            if logging:
                time_xrs += [mtg_to_dataset(g, variables=log_outputs, time=step + 1)]
                if len(time_xrs) >= max_time_steps_for_memory:
                    interstitial_dataset = xr.concat(time_xrs, dim="t")
                    interstitial_dataset.to_netcdf(output_path + f'/t={step + 1}.nc')
                    del interstitial_dataset
                    del time_xrs
                    time_xrs = []

            print(step)

        if logging:
            if len(time_xrs) > 0:
                interstitial_dataset = xr.concat(time_xrs, dim="t")
                interstitial_dataset.to_netcdf(output_path + f'/tf.nc')
                del interstitial_dataset
                del time_xrs

            # SAVING and merging
            # NOTE : merging is slower but way less space is needed
            time_step_files = [output_path + '/' + name for name in os.listdir(output_path)]
            time_step_datasets = [xr.load_dataset(file) for file in time_step_files]
            time_dataset = xr.concat(time_step_datasets, dim="t").sortby(["t", "vid"])
            # Note : you end up with nan when segment has not emerged

            time_dataset = time_dataset.assign_coords(coords=scenario).expand_dims(
                dim=dict(zip(list(scenario.keys()), [1 for k in scenario])))
            time_dataset.to_netcdf(output_path + '/merged.nc')
            del time_step_datasets
            del time_dataset
            for file in os.listdir(output_path):
                if '.nc' in file and file != "merged.nc":
                    os.remove(output_path + '/' + file)

        # Printing the time required by the simulation
        tmp = (time() - t_start) / 60.
        print("Simulation took %4.3f minutes!" % tmp)
        # Sound when simulation ends
        import winsound
        winsound.Beep(100, 500)

    else:
        # TODO
        print("You should save and reimport last step mtg in the output directory instead")
        root_bridges = Model(time_step=time_step)
        g = root_bridges.g

    time_dataset = xr.load_dataset(output_path + '/merged.nc')

    # Launching outputs analyses
    launch_analysis(dataset=time_dataset, mtg=g, output_dir=output_path,
                    global_state_extracts=global_state_extracts, global_flow_extracts=global_flow_extracts,
                    state_extracts=state_extracts, flow_extracts=flow_extracts,
                    global_sensitivity=False, global_plots=False, plot_architecture=True, STM_clustering=True)

    # TODO : Compare Cluster proportion of total flow or total pool at step t.
    # TODO : check units of states to verify before cumulating, te verify we talk about concentrations or quantity
