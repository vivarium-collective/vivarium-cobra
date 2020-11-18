# import registry
from vivarium.core.registry import (
    updater_registry,
    process_registry,
)

# import processes
from vivarium_cobra.processes.volume import Volume

# register processes
process_registry.register(Volume.name, Volume)


# register updaters
