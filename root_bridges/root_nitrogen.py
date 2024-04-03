from root_cynaps.root_nitrogen import RootNitrogenModel

from dataclasses import dataclass
from metafspm.component_factory import *
from metafspm.component import declare


@dataclass
class RootNitrogenModelCoupled(RootNitrogenModel):
    """

    NOTE : double names in methods are forbiden as they will be overwritten in Choregrapher's resolution
    However, it is a behavios of interest here when inheriting the class to edit it.
    """
    # INPUTS
    # FROM GROWTH MODEL
    amino_acids_consumption_by_growth: float = declare(default=0., unit="mol.s-1", unit_comment="", description="amino_acids consumption rate by growth processes", 
                                                 min_value="", max_value="", value_comment="", references="", DOI="",
                                                  variable_type="input", by="model_growth", state_variable_type="", edit_by="user")


    def __init__(self, g, time_step, **scenario):
        """Pass to inherited init, necessary with data classes"""
        super().__init__(g, time_step, **scenario)

    @state
    def _AA(self, AA, struct_mass, diffusion_AA_phloem, import_AA, diffusion_AA_soil, export_AA, AA_synthesis, 
            storage_synthesis, storage_catabolism, AA_catabolism, struct_mass_produced, amino_acids_consumption_by_growth):
        """
        EDIT : replaced structural nitrogen synthesis by rhizodep input in balance
        """
        if struct_mass > 0:
            return AA + (self.time_step / struct_mass) * (
                    diffusion_AA_phloem
                    + import_AA
                    - diffusion_AA_soil
                    - export_AA
                    + AA_synthesis
                    - storage_synthesis * self.r_AA_stor
                    + storage_catabolism / self.r_AA_stor
                    - AA_catabolism
                    - amino_acids_consumption_by_growth
            ) 
        else:
            return 0