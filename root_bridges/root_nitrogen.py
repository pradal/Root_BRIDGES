from root_cynaps.root_nitrogen import RootNitrogenModel

from dataclasses import dataclass
from metafspm.component_factory import *
from metafspm.component import declare

from numpy import pi, sqrt

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

    # STATE VARIABLES
    
    # PARAMETERS

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
            ) - struct_mass_produced * 0.2 / 146
        # glutamine 5 C -> 60g.mol-1 2N -> 28 g.mol-1 : C:N = 2.1
        # Sachant C:N struct environ de 10 = (Chex + CAA)/NAA Chex = 10*28 - 60 = 220 g Chex.
        # Sachang qu'un hexose contient 12*6=72 gC.mol-1 hex, c'est donc environ 3 hexoses pour 1 AA qui seraient consommés.
        # La proportion d'AA consommée par g de struct mass est donc de 1*146/(3*180 + 1*146) = 0.2 (180 g.mol-1 pour le glucose)

        else:
            return 0