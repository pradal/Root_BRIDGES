from rhizodep.soil_model import RhizoInputsSoilModel


class SoilModel(RhizoInputsSoilModel):
    """
    Empty doc
    """
    def __init__(self, g, time_step, **scenario):
        """Pass to inherited init, necessary with data classes"""
        super().__init__(g, time_step, **scenario)