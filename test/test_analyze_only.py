import os
from analyze.analyze import analyze_data
from log.logging import Logger


def test_anlyze_only():

    analyze_data(scenarios=["Reference"], 
                 outputs_dirpath="outputs/Reference", 
                 target_properties=None, **Logger.light_log)


test_anlyze_only()
