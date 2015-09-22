%{
connectivity.CellTestedPair (manual) # table of tested cell pair. Untested cell pair should not be included in the table
->connectivity.Slice
pair_id             : int            # id of a tested cell pair, directional
-----
connected           : tinyint        # 1 is connected, 0 is unconnected, 2 is EC coupling
%}

classdef CellTestedPair < dj.Relvar
end