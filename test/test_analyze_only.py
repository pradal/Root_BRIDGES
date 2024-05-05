from analyze.analyze import analyze_data


def test_anlyze_only():
    analyze_data(outputs_dirpath="outputs\T1",
                 on_sums=False,
                 on_performance=False,
                 animate_raw_logs=True,
                 target_properties=[]
                 )


test_anlyze_only()
