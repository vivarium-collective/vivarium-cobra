# import registry
from vivarium.core.registry import (
    updater_registry,
    process_registry,
)

# import processes
from vivarium_cobra.processes.volume import Volume
from vivarium_cobra.processes.local_field import LocalField
from vivarium_cobra.processes.homogeneous_environment import HomogeneousEnvironment

# register processes
process_registry.register(Volume.name, Volume)
process_registry.register(LocalField.name, LocalField)
process_registry.register(HomogeneousEnvironment.name, HomogeneousEnvironment)

# register updaters
