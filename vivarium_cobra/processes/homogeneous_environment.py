import numpy as np
from scipy import constants

from vivarium.library.units import units
from vivarium.core.process import Deriver

AVOGADRO = constants.N_A * 1 / units.mol


class HomogeneousEnvironment(Deriver):
    '''A non-spatial environment with volume'''

    name = 'homogeneous_environment'
    defaults = {
        'volume': 1e-12 * units.L,
        'concentrations': {},
    }

    def __init__(self, parameters=None):
        super(HomogeneousEnvironment, self).__init__(parameters)
        volume = parameters.get('volume', self.defaults['volume'])
        self.mmol_to_counts = (AVOGADRO.to('1/mmol') * volume).to('L/mmol')

    def ports_schema(self):
        return {
            'exchange': {
                '*': {
                    '_value': 0,
                },
            },
            'fields': {
                '*': {
                    '_default': 1.0,
                },
            },
        }



    def next_update(self, timestep, states):
        fields = states['fields']

        update = {
            'exchange': {
                mol_id: {
                    '_updater': 'set',
                    '_value': field[0][0],
                }
                for mol_id, field in fields.items()
            },
        }

        return update
