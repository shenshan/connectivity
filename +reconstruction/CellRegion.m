%{
reconstruction.CellRegion (lookup) # table for cell region labels
cell_region_id  : tinyint                # cell region identifier
---
cell_region_name            : char(8)                       # name of the cell region
%}


classdef CellRegion < dj.Relvar
end