'''
==================
Metabolism Process
==================

This module defines a :term:`process class` for modeling a cell's
metabolic processes with flux balance analysis (FBA). The cobrapy
FBA library is used for solving the problems. This supports metabolic
models from the `BiGG model database <http://bigg.ucsd.edu>`_,
and other configurations that can be passed to :py:class:`Metabolism`
to create models of metabolism.
'''

import os
import argparse

import numpy as np

from vivarium.core.process import Process
from vivarium.core.composition import (
    process_in_experiment,
    simulate_process,
    PROCESS_OUT_DIR,
)
from vivarium.plots.simulation_output import plot_simulation_output
from vivarium.library.units import units
from vivarium.library.dict_utils import tuplify_port_dicts

from vivarium_cobra.library.cobra_wrapper import FBA, AVOGADRO
from vivarium_cobra.library.regulation_logic import build_rule
from vivarium_cobra.processes.configurations import get_iAF1260b_config, get_toy_configuration

NAME = 'dynamic_fba'


def get_fg_from_counts(counts_dict, mw):
    composition_mass = sum([
        coeff / AVOGADRO * mw.get(mol_id, 0.0) * (units.g / units.mol)
        for mol_id, coeff in counts_dict.items()])  # g
    return composition_mass.to('fg')

fba_parameters = [
    'model_path',
    'default_tolerance',
    'tolerance',
    'stoichiometry',
    'reversible',
    'external_molecules',
    'objective',
    'flux_bounds',
    'molecular_weights',
    'exchange_bounds',
    'default_upper_bound',
    'target_added_mass',
]


class COBRA_FBA(Process):
    """A COBRA (COnstraint-Based Reconstruction and Analysis) process for metabolism

    This process runs flux balance analysis (FBA) with COBRApy.
    The FBA problem is defined using the provided configuration parameters,
    and can also load BiGG models.

    To see how to configure the process manually, look at the source
    code for :py:func:`test_toy_metabolism`.

    :term:`Ports`:

    * **external**: Holds the state of molecules external to the FBA
      reactions.
    * **internal_counts**: Holds the copy number of molecules internal
      to the FBA problem.
    * **exchange**: The accumulated counts of exchanged molecules with
      the external state, which can be used by an environmental compartment
      and updated for cell uptake and secretion.
    * **reactions**: Holds the names of the modeled metabolic reactions.
    * **flux_bounds**: Bounds on fluxes of the FBA problem.
    * **global**: Should be linked to the ``global`` :term:`store`.

    Args:
        initial_parameters (dict): Configures the process with the
            following keys/values:

            * **model_path** (:py:class:`str`): path to a BiGG model
            * **stoichiometry** (:py:class:`dict`): a map from reaction
              ID to that reaction's stoichiometry dictionary, e.g.
              ``{reaction_id: stoichiometry_dict}``
            * **objective** (:py:class:`dict`): the stoichiometry dict
              to be optimized
            * **external_molecules** (:py:class:`list`): the external
              molecules
            * **reversible_reactions** (:py:class:`list`)

    """

    name = NAME
    defaults = {
        'model_path': '',
        'stoichiometry': {},
        'reversible': [],
        'external_molecules': [],
        'objective': {},
        'flux_bounds': {},
        'molecular_weights': {},
        'exchange_bounds': {},
        'default_upper_bound': 0.0,
        'target_added_mass': 4.9e-7,  # approximates a doubling time of 2520 sec (42 min) in iAF1260b
        'regulation': {},
        'initial_state': {},
        'exchange_threshold': 1e-4,
        'initial_mass': 1000 * units.fg,
        'default_mmol_to_counts': 547467.333 * units.L / units.mmol,  # derived from 1000 fg mass
        'no_negative_objective': True,
        'time_step': 10,
    }

    def __init__(self, parameters=None):
        super().__init__(parameters)

        # initialize COBRA FBA
        self.fba = FBA({parameter: value
                        for parameter, value in self.parameters.items()
                        if parameter in fba_parameters})
        self.reaction_ids = self.fba.reaction_ids()

        # make the regulation functions
        regulation_logic = self.parameters['regulation']
        self.regulation = {
            reaction: build_rule(logic)
            for reaction, logic in regulation_logic.items()}

        # get internal molecules from fba objective
        self.objective_composition = self.fba.get_objective_composition()

        # configure initial mass with initial_state
        self.initial_mass = self.parameters['initial_mass']
        self.initial_external = {}
        self.initial_state()

    def initial_state(self, config=None):
        if config is None:
            config = {}
        initial_mass = config.get('initial_mass', self.parameters['initial_mass'])
        scale_concentration = config.get('scale_concentration', 1)
        override_initial = config.get('override_initial', {})

        # determine initial internal state from initial_mass
        mw = self.fba.molecular_weights
        composition = {
            mol_id: (-coeff if coeff < 0 else 0)
            for mol_id, coeff in self.objective_composition.items()}
        composition_mass = get_fg_from_counts(composition, mw)
        scaling_factor = (initial_mass / composition_mass).magnitude
        internal_state = {mol_id: int(coeff * scaling_factor)
            for mol_id, coeff in composition.items()}
        self.initial_mass = get_fg_from_counts(internal_state, mw)

        # Get external state from minimal_external fba solution
        self.initial_external = {
            state_id: 0.0
            for state_id in self.fba.external_molecules}
        self.initial_external.update(self.fba.minimal_external())

        # solve the fba problem to get flux_bounds
        exchange_fluxes = self.fba.read_exchange_fluxes()
        internal_fluxes = self.fba.read_internal_fluxes()
        flux_bounds = {
            mol: -val
            for mol, val in exchange_fluxes.items()}
        flux_bounds.update(internal_fluxes)

        # adjust external_state with config
        self.initial_external = {
            mol_id: conc * scale_concentration
            for mol_id, conc in self.initial_external.items()}
        for mol_id, conc in override_initial.items():
            self.initial_external[mol_id] = conc

        return {
            'external': self.initial_external,
            'internal_counts': internal_state}

    def ports_schema(self):
        ports = [
            'internal_counts',
            'external',
            'exchanges',
            'reactions',
            'flux_bounds',
            'global']

        schema = {port: {} for port in ports}

        # internal counts
        for state in list(self.objective_composition.keys()):
            schema['internal_counts'][state] = {
                '_divider': 'split',
                '_default': 0.0,
                '_emit': True,
                '_properties': {
                    'mw': self.fba.molecular_weights[state] * units.g / units.mol}}

        # external and exchange
        for state in self.fba.external_molecules:
            schema['external'][state] = {
                '_default': self.initial_external[state],
                '_emit': True}
            schema['exchanges'][state] = {
                '_default': 0}

        # reactions
        for state in self.reaction_ids:
            schema['reactions'][state] = {
                '_default': 0.0,
                '_updater': 'set'}

        # flux bounds
        schema['flux_bounds'] = {
            '*': {
                '_default': self.parameters['default_upper_bound'],
                '_emit': True}}

        # globals
        schema['global'] = {
            'mass': {
                '_default': self.initial_mass,
                '_emit': True},
            'mmol_to_counts': {
                '_default': self.parameters['default_mmol_to_counts'],
            }}

        return schema


    def next_update(self, timestep, states):

        # get the state
        external_state = states['external']
        flux_bounds = states['flux_bounds']  # mmol/L/s
        mmol_to_counts = states['global']['mmol_to_counts'].to('L/mmol').magnitude

        # get constraints
        ## exchange_constraints based on external availability
        exchange_constraints = {mol_id: 0.0
            for mol_id, conc in external_state.items()
                                if conc <= self.parameters['exchange_threshold']}

        ## state of regulated reactions (True/False)
        flattened_states = tuplify_port_dicts(states)
        regulation_state = {}
        for reaction_id, reg_logic in self.regulation.items():
            regulation_state[reaction_id] = reg_logic(flattened_states)

        # apply constraints
        ## exchange constraints
        self.fba.set_exchange_bounds(exchange_constraints)

        ## constraints from flux_bounds
        if flux_bounds:
            # only pass in reaction_ids that exist in the fba model
            constrained_reaction_bounds = {
                reaction_id: constraint
                for reaction_id, constraint in flux_bounds.items()
                if reaction_id in self.reaction_ids}
            self.fba.constrain_flux(constrained_reaction_bounds)

        ## turn reactions on/off based on regulation
        self.fba.regulate_flux(regulation_state)

        # solve the fba problem
        objective_exchange = self.fba.optimize() * timestep  # mmol/L/s
        exchange_reactions = self.fba.read_exchange_reactions()
        exchange_fluxes = self.fba.read_exchange_fluxes()  # mmol/L/s
        internal_fluxes = self.fba.read_internal_fluxes()  # mmol/L/s

        # if no valid solution, return empty update
        if (self.parameters['no_negative_objective'] and objective_exchange < 0.0) \
                or np.isnan(objective_exchange):
            return {}

        # convert results
        ## time step dependence on fluxes
        exchange_fluxes.update((mol_id, flux * timestep)
                               for mol_id, flux in exchange_fluxes.items())
        internal_fluxes.update((mol_id, flux * timestep)
                               for mol_id, flux in internal_fluxes.items())

        ## update internal counts from objective flux
        ## calculate added mass from the objective molecules' molecular weights
        objective_count = objective_exchange * mmol_to_counts

        internal_state_update = {}
        for reaction_id, coeff1 in self.fba.objective.items():
            for mol_id, coeff2 in self.fba.stoichiometry[reaction_id].items():
                if coeff2 < 0:  # pull out molecule if it is used to make biomass (negative coefficient)
                    added_count = int(-coeff1 * coeff2 * objective_count)
                    internal_state_update[mol_id] = added_count

        ## convert exchange fluxes to counts
        exchanges_updates = {
            mol_id: int(flux * mmol_to_counts)
            for mol_id, flux in exchange_fluxes.items()}

        all_fluxes = {}
        all_fluxes.update(internal_fluxes)
        all_fluxes.update(exchange_reactions)

        update = {
            'exchanges': exchanges_updates,
            'internal_counts': internal_state_update,
            'reactions': all_fluxes}

        return update


# tests
def test_toy_metabolism(
        total_time=15
):
    regulation_logic = {
        'R4': 'if (external, O2) > 0.1 and not (external, F) < 0.1'}
    toy_config = get_toy_configuration()
    toy_config['regulation'] = regulation_logic
    toy_config['target_added_mass'] = None
    toy_metabolism = COBRA_FBA(toy_config)

    # simulate toy model
    interval = int(total_time/3)
    timeline = [
        (interval, {('external', 'A'): 1}),
        (2*interval, {('external', 'F'): 0}),
        (total_time, {})]

    settings = {
        'initial_state': toy_metabolism.initial_state(),
        # 'total_time': total_time,
        'timeline': {
            'timeline': timeline,
            'time_step': 10}
    }
    return simulate_process(toy_metabolism, settings)


def run_bigg(
        total_time=100,
        volume=1e-5,
):
    config = get_iAF1260b_config()
    metabolism = COBRA_FBA(config)
    initial_config = {}
    initial_state = metabolism.initial_state(
        config=initial_config)

    # run simulation
    settings = {
        'initial_state': initial_state}
    experiment = process_in_experiment(metabolism, settings)

    # run simulation
    experiment.update(total_time)
    experiment.end()

    # return data from emitter
    return experiment.emitter.get_timeseries()


def print_growth(global_timeseries):
    volume_ts = global_timeseries[('volume', 'femtoliter')]
    mass_ts = global_timeseries[('mass', 'femtogram')]
    print('volume growth: {}'.format(volume_ts[-1] / volume_ts[0]))
    print('mass growth: {}'.format(mass_ts[-1] / mass_ts[0]))


def main():
    out_dir = os.path.join(PROCESS_OUT_DIR, NAME)
    os.makedirs(out_dir, exist_ok=True)

    parser = argparse.ArgumentParser(description='metabolism process')
    parser.add_argument('--bigg', '-b', action='store_true', default=False, )
    parser.add_argument('--toy', '-t', action='store_true', default=False, )
    args = parser.parse_args()

    if args.toy:
        timeseries = test_toy_metabolism(
            total_time=2500)

        print_growth(timeseries['global'])

        plot_settings = {}
        plot_simulation_output(timeseries, plot_settings, out_dir, 'toy_metabolism')

    else:
        timeseries = run_bigg(
            total_time=2500)

        # save_timeseries(timeseries, out_dir)  # TODO -- make a test with timeseries reference
        print_growth(timeseries['global'])

        # plot
        plot_settings = {
            'max_rows': 30,
            'remove_zeros': True,
            'skip_ports': ['exchange', 'reactions']}
        plot_simulation_output(timeseries, plot_settings, out_dir, 'BiGG_simulation')


if __name__ == '__main__':
    main()
