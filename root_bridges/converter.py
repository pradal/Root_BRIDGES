class Assembly:
    def __init__(self):
        pass

    @property
    def link(self):
        # As variables to sum have already been added to self during link_mtg call...
        dict_list = [d for d in dir(self) if ("__" not in d and d != "link")]
        output_dict = {}
        converted_dicts = []
        for d in dict_list:
            # For all extracted dicts
            working_dict = getattr(self, d)
            # We list all dicts after conversion factors have been applied to prepare them for sum
            converted_dicts += [dict(zip(working_dict["values"].keys(), [v*working_dict["conversion"]
                                                                         for v in working_dict["values"].values()]))]
        # for each "vid"
        for key in converted_dicts[0].keys():
            total = 0
            # for each property to sum
            for d in converted_dicts:
                # we append the total
                total += d[key]
            output_dict.update({key: total})
        # retrieve the actualized value on each property call in the receiver class
        return output_dict


def link_mtg(receiver, applier, category, translator={}, same_names=True):
    """
    Description : linker function that will enable properties sharing through MTG.

    Parameters :
    :param receiver: (class) model class whose inputs should be provided with the applier class.
    :param applier: (class) model class whose properties are used to provide inputs to the receiver class.
    :param category: (sting) word to specify which inputs are to be considered in the receiver model class.
    :param translator: (dict) translation dict used when receiver and applier properties do not have the same names.
    :param same_names: (bool) boolean value to be used if a model was developped by another team with different names.

    Note :  The whole property is transfered, so if only the collar value of a spatial property is needed,
    it will be accessed through the first vertice with the [1] indice. Not spatialized properties like xylem pressure or
    single point properties like collar flows are only stored in the indice [1] vertice.
    """
    if same_names:
        for link in getattr(receiver, "inputs")[category]:
            setattr(receiver, link, getattr(applier, link))
    else:
        for link in getattr(receiver, "inputs")[category]:
            # We create an instance of the above class with the name of the target property we build from applier inputs
            setattr(receiver, link, Assembly())
            # For each pointed variables
            for prop in translator[link]:
                # We link them to the nested instantiated class
                setattr(getattr(receiver, link), prop, {"values": getattr(applier, prop), "conversion": translator[link][prop]})
                # Then they will be converted and summed on call of the Assembly.link property


# Example dict to align different names, conversion factors, and properties to aggregate
shared_states = {
    "external_property": {
        "pool_3": 1,  # comment Here unit -> unit (hypotheses)
        "pool_4": 1   # comment Here unit -> unit (hypotheses)
    },
    "external_property_2": {
        "pool_3": 10
    }
}
