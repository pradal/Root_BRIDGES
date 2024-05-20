import os
# Utility packages
from analyze.analyze import analyze_data


if __name__ == '__main__':

    scenario_name = "Reference_Fischer"

    analyze_data(outputs_dirpath=os.path.join("outputs", scenario_name),
                     on_sums=True,
                     on_performance=True,
                     animate_raw_logs=True,
                     target_properties=None
                     )
    # In the end put the system to sleep, Windows only
    #os.system("rundll32.exe powrprof.dll,SetSuspendState 0,1,0")
