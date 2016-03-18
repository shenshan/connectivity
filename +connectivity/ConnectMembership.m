%{
connectivity.ConnectMembership (manual) # my newest table
-> connectivity.CellTestedPair
-> connectivity.Cell
-----
role      : enum('from', 'to', 'EC')        # specify the role of the cell, from means the cell is stimulated, EC represents EC coupling
%}

classdef ConnectMembership < dj.Relvar

end