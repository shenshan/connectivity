%{
reconstruction.DensityParameters (lookup) # 
density_param_id: tinyint                # identifier of the parameter set
---
bin_width                   : float                         # bin width for the 3d histogram in mu
raster_resolution           : float                         # resolution of the rasterization in mu
cell_offset                 : float                         # assumed offset between cells into the slice
multiplication_factor       : int                           # multiply the number of points by this factor when spin shuffling
cut_compensation            : tinyint unsigned              # try to account for the fact that parts of the neuron might be missing
%}


classdef DensityParameters < dj.Relvar
end