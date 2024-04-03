from rhizodep.root_growth import RootGrowthModel

from dataclasses import dataclass
from metafspm.component_factory import *
from metafspm.component import declare

@dataclass
class RootGrowthModelCoupled(RootGrowthModel):
    """

    NOTE : double names in methods are forbiden as they will be overwritten in Choregrapher's resolution
    However, it is a behavios of interest here when inheriting the class to edit it.
    """
    # INPUTS
    
    # STATE VARIABLES
    
    # PARAMETERS
    


    def __init__(self, g, time_step, **scenario):
        """Pass to inherited init, necessary with data classes"""
        super().__init__(g, time_step, **scenario)


    # Function for calculating root elongation:
    def elongated_length(self, element, initial_length: float, radius: float, C_hexose_root: float, elongation_time_in_seconds: float):
        """
        This function computes a new length (m) based on the elongation process described by ArchiSimple and regulated by
        the available concentration of hexose.
        It has been modified in Root-BRIDGES to introduce a regulation by 
        :param initial_length: the initial length (m)
        :param radius: radius (m)
        :param C_hexose_root: the concentration of hexose available for elongation (mol of hexose per gram of strctural mass)
        :param elongation_time_in_seconds: the period of elongation (s)
        :return: the new elongated length
        """

        # If we keep the classical ArchiSimple rule:
        if self.ArchiSimple:
            # Then the elongation is calculated following the rules of Pages et al. (2014):
            elongation = self.EL * 2. * radius * elongation_time_in_seconds
        else:
            # Otherwise, we additionally consider a limitation of the elongation according to the local concentration of hexose,
            # based on a Michaelis-Menten formalism:
            if C_hexose_root > 0.:
                elongation = self.EL * 2. * radius / (
                    ((1 + self.Km_elongation) / C_hexose_root) * ((1 + self.Km_elongation) / element.AA)
                        ) * elongation_time_in_seconds
            else:
                elongation = 0.

        # We calculate the new potential length corresponding to this elongation:
        new_length = initial_length + elongation
        if new_length < initial_length:
            print("!!! ERROR: There is a problem of elongation, with the initial length", initial_length,
                  " and the radius", radius, "and the elongation time", elongation_time_in_seconds)
        return new_length
    

    # Function for calculating the amount of C to be used in neighbouring elements for sustaining root elongation:
    def calculating_supply_for_elongation(self, element):
        """
        This function computes the list of root elements that can supply C as hexose for sustaining the elongation
        of a given element, as well as their structural mass and their amount of available hexose.
        :param element: the element for which we calculate the possible supply of C for its elongation
        :return: three lists containing the indices of elements, their hexose amount (mol of hexose) and their structural mass (g).
        """
        # TODO FOR TRISTAN: Consider using a similar approach for sustaining the need for N when a root apex should elongate.
        n = element

        # We initialize each amount of hexose available for growth:
        n.hexose_possibly_required_for_elongation = 0.
        n.amino_acids_possibly_required_for_elongation = 0.
        n.struct_mass_contributing_to_elongation = 0.

        # We initialize empty lists:
        list_of_elongation_supporting_elements = []
        list_of_elongation_supporting_elements_hexose = []
        list_of_elongation_supporting_elements_amino_acids = []
        list_of_elongation_supporting_elements_mass = []

        # We then calculate the length of an apical zone of a fixed length which can provide the amount of hexose required for growth:
        growing_zone_length = self.growing_zone_factor * n.radius
        # We calculate the corresponding volume to which this length should correspond based on the diameter of this apex:
        supplying_volume = growing_zone_length * n.radius ** 2 * pi

        # We start counting the hexose at the apex:
        index = n.index()
        current_element = n

        # We initialize a temporary variable that will be used as a counter:
        remaining_volume = supplying_volume

        # As long the remaining volume is not zero:
        while remaining_volume > 0:

            # If the volume of the current element is lower than the remaining volume:
            if remaining_volume > current_element.volume:
                # We make sure to include in the list of supplying elements only elements with a positive length
                # (e.g. NOT the elements of length 0 that support seminal or adventitious roots):
                if current_element.length > 0.:
                    # We add to the amount of hexose available all the hexose in the current element
                    # (EXCLUDING sugars in the living root hairs):
                    # TODO: Should the C from root hairs be used for helping roots to grow?
                    hexose_contribution = current_element.C_hexose_root * current_element.struct_mass
                    amino_acids_contribution = current_element.AA * current_element.struct_mass
                    n.hexose_possibly_required_for_elongation += hexose_contribution
                    n.amino_acids_possibly_required_for_elongation += amino_acids_contribution
                    n.struct_mass_contributing_to_elongation += current_element.struct_mass
                    # We record the index of the contributing element:
                    list_of_elongation_supporting_elements.append(index)
                    # We record the amount of hexose that the current element can provide:
                    list_of_elongation_supporting_elements_hexose.append(hexose_contribution)
                    # We record the amount of amino acids that the current element can provide:
                    list_of_elongation_supporting_elements_amino_acids.append(amino_acids_contribution)
                    # We record the structural mass from which the current element contributes:
                    list_of_elongation_supporting_elements_mass.append(current_element.struct_mass)
                    # We subtract the volume of the current element to the remaining volume:
                    remaining_volume = remaining_volume - current_element.volume

                # And we try to move the index to the segment preceding the current element:
                index_attempt = self.g.Father(index, EdgeType='<')
                # If there is no father element on this axis:
                if index_attempt is None:
                    # Then we try to move to the mother root, if any:
                    index_attempt = self.g.Father(index, EdgeType='+')
                    # If there is no such root:
                    if index_attempt is None:
                        # Then we exit the loop here:
                        break
                # We set the new index:
                index = index_attempt
                # We define the new element to consider according to the new index:
                current_element = self.g.node(index)
            # Otherwise, this is the last preceding element to consider:
            else:
                # We finally add to the amount of hexose available for elongation a part of the hexose of the current element:
                hexose_contribution = current_element.C_hexose_root * current_element.struct_mass \
                                      * remaining_volume / current_element.volume
                amino_acids_contribution = current_element.AA * current_element.struct_mass \
                                      * remaining_volume / current_element.volume
                n.hexose_possibly_required_for_elongation += hexose_contribution
                n.amino_acids_possibly_required_for_elongation += amino_acids_contribution
                n.struct_mass_contributing_to_elongation += current_element.struct_mass \
                                                            * remaining_volume / current_element.volume
                # We record the index of the contributing element:
                list_of_elongation_supporting_elements.append(index)
                # We record the amount of hexose that the current element can provide:
                list_of_elongation_supporting_elements_hexose.append(hexose_contribution)
                # We record the amount of amino acids that the current element can provide:
                list_of_elongation_supporting_elements_amino_acids.append(amino_acids_contribution)
                # We record the structural mass from which the current element contributes:
                list_of_elongation_supporting_elements_mass.append(
                    current_element.struct_mass * remaining_volume / current_element.volume)
                # And the remaining volume to consider is set to 0:
                remaining_volume = 0.
                # And we exit the loop here:
                break

        # We record the average concentration in hexose of the whole zone of hexose supply contributing to elongation:
        if n.struct_mass_contributing_to_elongation > 0.:
            n.growing_zone_C_hexose_root = n.hexose_possibly_required_for_elongation / n.struct_mass_contributing_to_elongation
        else:
            print("!!! ERROR: the mass contributing to elongation in element", n.index(), "of type", n.type, "is",
                  n.struct_mass_contributing_to_elongation,
                  "g, and its structural mass is", n.struct_mass, "g!")
            n.growing_zone_C_hexose_root = 0.

        n.list_of_elongation_supporting_elements = list_of_elongation_supporting_elements
        n.list_of_elongation_supporting_elements_hexose = list_of_elongation_supporting_elements_hexose

        n.list_of_elongation_supporting_elements_mass = list_of_elongation_supporting_elements_mass

        return list_of_elongation_supporting_elements, \
            list_of_elongation_supporting_elements_hexose, \
            list_of_elongation_supporting_elements_mass