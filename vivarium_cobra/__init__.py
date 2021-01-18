# import registry
from vivarium.core.registry import process_registry

# import processes
from vivarium_cobra.processes.volume import Volume
from vivarium_cobra.processes.local_field import LocalField

# register processes
process_registry.register(Volume.name, Volume)
process_registry.register(LocalField.name, LocalField)
