

from vivarium.core.process import Composite
from vivarium.processes.tree_mass import TreeMass
from vivarium_cobra.processes.volume import Volume
from vivarium_cobra.processes.dynamic_fba import DynamicFBA
from vivarium_cobra.processes.configurations import get_iAF1260b_config


class CobraComposite(Composite):

    # set cobra constrained reactions
    cobra_config = get_iAF1260b_config()

    defaults = {
        'cobra': cobra_config,
        'flux_deriver': {},
    }

    def generate_processes(self, config):

        cobra_process = DynamicFBA(config['cobra'])

        return {
            'cobra': cobra_process,
            'mass_deriver': TreeMass(),
            'volume_deriver': Volume(),
        }

    def generate_topology(self, config):
        globals_path = ('globals',)

        return {
            'cobra': {
                'internal_counts': ('internal_counts',),
                'external': ('external',),
                'exchanges': ('exchanges',),
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


def main():

    import ipdb; ipdb.set_trace()




if __name__ == '__main__':
    main()
