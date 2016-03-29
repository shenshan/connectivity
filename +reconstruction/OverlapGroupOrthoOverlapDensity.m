%{
reconstruction.OverlapGroupOrthoOverlapDensity (computed) # overlap densities for a particular cell type pair and distance
-> reconstruction.OverlapGroup
-> connectivity.Distance
---
d                           : longblob                      # distance to the connecting line between two neurons
p                           : longblob                      # overlap density along that line
%}


classdef OverlapGroupOrthoOverlapDensity < dj.Relvar
	methods

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			self.insert(key)
		end
	end

end