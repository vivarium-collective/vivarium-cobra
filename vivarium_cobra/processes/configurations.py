import os


def get_package_path():
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


def get_e_coli_core_config():
    """Get an *E. coli* core metabolism model

    The model is the `e_coli_core model from BiGG
    <http://bigg.ucsd.edu/models/e_coli_core>`_.

    Returns:
        A configuration for the model that can be passed to the
        :py:class:`Metabolism` constructor.
    """
    package_path = get_package_path()
    metabolism_file = os.path.join(package_path, 'bigg_models', 'e_coli_core.json')
    return {'model_path': metabolism_file}


def get_iAF1260b_config():
    """Get the metabolism config for the iAF1260b BiGG model

    The metabolism model is the `iAF1260b model from BiGG
    <http://bigg.ucsd.edu/models/iAF1260b>`_.

    Returns:
        A configuration for the model that can be passed to the
        :py:class:`Metabolism` constructor.
    """
    package_path = get_package_path()
    metabolism_file = os.path.join(package_path, 'bigg_models', 'iAF1260b.json')
    return {'model_path': metabolism_file}


def get_toy_configuration():
    stoichiometry = {
        'R1': {'A': -1, 'ATP': -1, 'B': 1},
        'R2a': {'B': -1, 'ATP': 2, 'NADH': 2, 'C': 1},
        'R2b': {'C': -1, 'ATP': -2, 'NADH': -2, 'B': 1},
        'R3': {'B': -1, 'F': 1},
        'R4': {'C': -1, 'G': 1},
        'R5': {'G': -1, 'C': 0.8, 'NADH': 2},
        'R6': {'C': -1, 'ATP': 2, 'D': 3},
        'R7': {'C': -1, 'NADH': -4, 'E': 3},
        'R8a': {'G': -1, 'ATP': -1, 'NADH': -2, 'H': 1},
        'R8b': {'G': 1, 'ATP': 1, 'NADH': 2, 'H': -1},
        'Rres': {'NADH': -1, 'O2': -1, 'ATP': 1},
        'v_biomass': {'C': -1, 'F': -1, 'H': -1, 'ATP': -10}}

    external_molecules = ['A', 'F', 'D', 'E', 'H', 'O2']

    objective = {'v_biomass': 1.0}

    reversible = ['R6', 'R7', 'Rres']

    default_reaction_bounds = 1000.0

    exchange_bounds = {
        'A': -0.02,
        'D': 0.01,
        'E': 0.01,
        'F': -0.005,
        'H': -0.005,
        'O2': -0.1}

    initial_state = {
        'external': {
            'A': 21.0,
            'F': 5.0,
            'D': 12.0,
            'E': 12.0,
            'H': 5.0,
            'O2': 100.0}}

    # molecular weight units are (units.g / units.mol)
    molecular_weights = {
        'A': 500.0,
        'B': 500.0,
        'C': 500.0,
        'D': 500.0,
        'E': 500.0,
        'F': 50000.0,
        'H': 1.00794,
        'O2': 31.9988,
        'ATP': 507.181,
        'NADH': 664.425}

    config = {
        'stoichiometry': stoichiometry,
        'reversible': reversible,
        'external_molecules': external_molecules,
        'objective': objective,
        'initial_state': initial_state,
        'exchange_bounds': exchange_bounds,
        'default_upper_bound': default_reaction_bounds,
        'molecular_weights': molecular_weights}

    return config
