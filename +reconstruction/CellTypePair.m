%{
reconstruction.CellTypePair (lookup) # table with cell type name pairs fetched from CellType
cell_type_from  : varchar(256)           # 
cell_type_to    : varchar(256)           # 
---
%}


classdef CellTypePair < dj.Relvar
end