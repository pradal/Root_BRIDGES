from rhizodep.soil_model import RhizoInputsSoilModel

from dataclasses import dataclass
from metafspm.component_factory import *
from metafspm.component import declare

@dataclass
class SoilModel(RhizoInputsSoilModel):
    """
    Empty doc
    """
    # INPUTS
    # FROM NITROGEN MODEL
    mineralN_uptake: float = declare(default=0., unit="mol.s-1", unit_comment="of nitrates", description="", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="", 
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    amino_acids_uptake: float = declare(default=0., unit="mol.s-1", unit_comment="of amino acids", description="", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    mineralN_diffusion_from_roots: float =  declare(default=0., unit="mol.s-1", unit_comment="of nitrates", description="", 
                                                    min_value="", max_value="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    amino_acids_diffusion_from_roots: float =  declare(default=0., unit="mol.s-1", unit_comment="of amino acids", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    mineralN_diffusion_from_xylem: float =  declare(default=0., unit="mol.s-1", unit_comment="of nitrates", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    amino_acids_diffusion_from_xylem: float =  declare(default=0., unit="mol.s-1", unit_comment="of amino_acids", 
                                                    min_value="", max_value="", description="", value_comment="", references="", DOI="",
                                                    variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")


    # STATE VARIABLES

    # PARAMETERS

    def __init__(self, g, time_step, **scenario):
        """Pass to inherited init, necessary with data classes"""
        super().__init__(g, time_step, **scenario)


    @state
    def _C_mineralN_soil(self, C_mineralN_soil, volume_soil,mineralN_diffusion_from_roots, mineralN_diffusion_from_xylem, mineralN_uptake):
        balance = C_mineralN_soil + (self.time_step_in_seconds / volume_soil) * (
            mineralN_diffusion_from_roots
            + mineralN_diffusion_from_xylem
            - mineralN_uptake
        )
        balance[balance < 0.] = 0.
        return balance

    @state
    def _C_amino_acids_soil(self, C_amino_acids_soil, volume_soil, amino_acids_diffusion_from_roots, amino_acids_diffusion_from_xylem, amino_acids_uptake):
        balance = C_amino_acids_soil + (self.time_step_in_seconds / volume_soil) * (
            amino_acids_diffusion_from_roots
            + amino_acids_diffusion_from_xylem
            - amino_acids_uptake
        )
        balance[balance < 0.] = 0.
        return balance