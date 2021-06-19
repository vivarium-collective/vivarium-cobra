import os

from vivarium.core.composer import Composer
from vivarium.core.composition import composer_in_experiment, COMPOSITE_OUT_DIR
from vivarium.library.units import units

# core processes
from vivarium.processes.tree_mass import TreeMass

# cobra processes
from vivarium_cobra.processes.cobra_fba import COBRA_FBA, print_growth
from vivarium_cobra.processes.volume import Volume
from vivarium_cobra.processes.local_field import LocalField
from vivarium_cobra.processes.configurations import get_iAF1260b_config

from vivarium.plots.simulation_output import plot_simulation_output



NAME = 'cobra_composite'

class CobraComposite(Composer):

    defaults = {
        'cobra': {},
        'flux_deriver': {},
        'fields_on': True,
        'field_deriver':   {
            'nonspatial': True,
            'bin_volume': 1e-6 * units.L}
    }

    def generate_processes(self, config):

        cobra_process = COBRA_FBA(config['cobra'])

        processes = {
            'cobra': cobra_process,
            'mass_deriver': TreeMass(),
            'volume_deriver': Volume(),
        }
        if config['fields_on']:
            initial_state_cobra = cobra_process.initial_state()
            fields_config = config['field_deriver']
            fields_config.update({
                'initial_external': initial_state_cobra['external']
            })
            processes.update({
                'field_deriver': LocalField(fields_config)
            })

        return processes

    def generate_topology(self, config):
        globals_path = ('global',)
        exchanges_path = ('exchanges',)
        external_path = ('external',)
        topology = {
            'cobra': {
                'internal_counts': ('internal_counts',),
                'external': external_path,
                'exchanges': exchanges_path,
                'reactions': ('reactions',),
                'flux_bounds': ('flux_bounds',),
                'global': globals_path,
            },
            'mass_deriver': {
                'global': globals_path,
            },
            'volume_deriver': {
                'global': globals_path,
            },
        }

        if config['fields_on']:
            topology.update({
                'field_deriver': {
                    'exchanges': exchanges_path,
                    'location': ('global', 'location',),
                    'fields': external_path,
                    'dimensions': ('dimensions',),
                    }
            })

        return topology

def test_cobra_composite(
        total_time=100,
):
    # initialize composite
    cobra_config = get_iAF1260b_config()
    config = {
        'cobra': cobra_config
    }
    composite = CobraComposite(config)

    # initial state
    initial_config = {}
    initial_state = composite.initial_state(
        config=initial_config)

    # make the experiment
    exp_settings = {
        'initial_state': initial_state,
    }
    experiment = composer_in_experiment(composite, exp_settings)

    # run
    experiment.update(total_time)

    # get the data
    return experiment.emitter.get_timeseries()


def plot_sim_output(timeseries, out_dir='out'):
    print_growth(timeseries['global'])

    # plot
    plot_settings = {
        'max_rows': 30,
        'remove_zeros': True,
        'skip_ports': ['exchange', 'reactions']}
    plot_simulation_output(
        timeseries,
        plot_settings,
        out_dir,
        'cobra_composite')


def main(total_time=100, out_dir='out'):
    timeseries = test_cobra_composite(
        total_time=total_time)
    plot_sim_output(timeseries, out_dir)

if __name__ == '__main__':
    out_dir = os.path.join(COMPOSITE_OUT_DIR, NAME)
    os.makedirs(out_dir, exist_ok=True)
    main(
        total_time=2500,
        out_dir=out_dir)
