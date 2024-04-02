from rhizodep.root_carbon import RootCarbonModel

from dataclasses import dataclass
from metafspm.component_factory import *
from metafspm.component import declare

@dataclass
class RootCarbonModelCoupled(RootCarbonModel):
    """

    NOTE : double names in methods are forbiden as they will be overwritten in Choregrapher's resolution
    However, it is a behavios of interest here when inheriting the class to edit it.
    """
    # INPUTS
    import_Nm: float =  declare(default=0., unit="mol.s-1", unit_comment="of nitrates", description="", 
                                min_value="", max_value="", value_comment="", references="", DOI="", 
                                variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    import_AA: float =  declare(default=0., unit="mol.s-1", unit_comment="of amino acids", description="", 
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    export_Nm: float =  declare(default=0., unit="mol.s-1", unit_comment="of nitrates", description="", 
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    export_AA: float =  declare(default=0., unit="mol.s-1", unit_comment="of amino acids", description="", 
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    AA_synthesis: float = declare(default=0., unit="mol.s-1", unit_comment="of amino acids", description="", 
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    AA_catabolism: float = declare(default=0., unit="mol.s-1", unit_comment="of amino acids", description="", 
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="input", by="model_nitrogen", state_variable_type="extensive", edit_by="user")
    
    # STATE VARIABLES
    N_metabolic_respiration: float = declare(default=0., unit="mol.s-1", unit_comment="of carbon", description="Respiration related to nitrogen metabolism", 
                                            min_value="", max_value="", value_comment="", references="", DOI="",
                                             variable_type="state_variable", by="model_carbon", state_variable_type="extensive", edit_by="user")

    # PARAMETERS
    r_hexose_AA: float = declare(default=4.5/6, unit="adim", unit_comment="mol of hexose per mol of amino acids in roots", description="stoechiometric ratio during amino acids synthesis for hexose consumption", 
                                min_value="", max_value="", value_comment="", references="we hypothesize from Yemm and Willis 1956 that synthetized soluble amino acids are mainly composed of glutamine and asparagine", DOI="",
                                variable_type="parameter", by="model_carbon", state_variable_type="", edit_by="user")
    r_Nm_AA: float =     declare(default=1.4, unit="adim", unit_comment="mol of N per mol of amino acids", description="concentration stoechiometric ratio between mineral nitrogen and amino acids in roots", 
                                min_value="", max_value="", value_comment="", references="", DOI="",
                                variable_type="parameter", by="model_nitrogen", state_variable_type="", edit_by="user")
    respi_costs_mineralN_import: float = declare(default=0.397, unit="adim", unit_comment="mol of C per mol of N", description="Respiratory of active imports in root.", 
                                min_value="", max_value="", value_comment="", references="Barillot et al., 2016", DOI="",
                                variable_type="parameter", by="model_carbon", state_variable_type="", edit_by="user")
    respi_costs_mineralN_reduction: float = declare(default=1.98, unit="adim", unit_comment="mol of C per mol of N", description="Respiratory of active imports in root.", 
                                min_value="", max_value="", value_comment="", references="Robinson 2001; Barillot et al., 2016", DOI="",
                                variable_type="parameter", by="model_carbon", state_variable_type="", edit_by="user")


    def __init__(self, g, time_step, **scenario):
        """Pass to inherited init, necessary with data classes"""
        super().__init__(g, time_step, **scenario)


    @rate
    def _N_metabolic_respiration(self, import_Nm, export_Nm, import_AA, export_AA, AA_synthesis):
        """
        Here we explicit a respiration cost associated to the exchange of N in the root (e.g. cost for uptake and xylem loading + anabolism costs of soluble nitrogen).        
        In line with Thornley and Cannell 2000; Robinson 2001 and Barillot et al. 2016, respiratory costs are considered as proportionnal to transport and synthesis costs.
        Moreover, such synthesis can have side effects of carbohydrates catabolism affecting the respiration through the release of labile hexoses (Yemm and Willis 1956).
        """
        print("respiring")
        transport_respiration = self.respi_costs_mineralN_import * (import_Nm + export_Nm + import_AA + export_AA)
        anabolism_respiration = self.respi_costs_mineralN_reduction * AA_synthesis *self.r_Nm_AA
        return transport_respiration + anabolism_respiration

    @potential
    @state
    def _C_hexose_root(self, C_hexose_root, struct_mass, living_root_hairs_struct_mass, hexose_exudation, hexose_uptake_from_soil,
                           mucilage_secretion, cells_release, maintenance_respiration,
                           hexose_consumption_by_growth, hexose_diffusion_from_phloem,
                           hexose_active_production_from_phloem, sucrose_loading_in_phloem,
                           hexose_mobilization_from_reserve, hexose_immobilization_as_reserve, deficit_hexose_root, 
                           AA_synthesis, AA_catabolism, N_metabolic_respiration):
        """
        Added the following flows to the balance :
        - Amino acid synthesis hexose consumption
        - Amino acid catabolism releasing hexose
        - Nitrogen metabolism related respiration costs
        """
        return C_hexose_root + (self.time_steps_in_seconds / (struct_mass + living_root_hairs_struct_mass)) * (
                - hexose_exudation
                + hexose_uptake_from_soil
                - mucilage_secretion
                - cells_release
                - maintenance_respiration / 6.
                - hexose_consumption_by_growth
                + hexose_diffusion_from_phloem
                + hexose_active_production_from_phloem
                - 2. * sucrose_loading_in_phloem
                + hexose_mobilization_from_reserve
                - hexose_immobilization_as_reserve
                - deficit_hexose_root
                - AA_synthesis * self.r_hexose_AA
                + AA_catabolism / self.r_hexose_AA
                - N_metabolic_respiration / 6.)