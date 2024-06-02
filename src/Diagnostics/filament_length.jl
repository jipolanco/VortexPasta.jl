# Note: the filament_length function used to be defined here in Diagnostics, but it was
# later moved to the Filaments module. We still keep the interface below mostly because it
# this should also be considered a diagnostic along with the energy and other stuff, so that
# it appears in the same docs page and can be extended via Diagnostics.filament_length.

using ..Filaments: filament_length
export filament_length
