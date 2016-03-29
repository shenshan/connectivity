%{
reconstruction.PaperConnectivity (lookup) # 
cell_type_from  : varchar(256)           # 
cell_type_to    : varchar(256)           # 
---
k                           : int                           # 
n                           : int                           # 
p                           : float                         # 
%}


classdef PaperConnectivity < dj.Relvar
end