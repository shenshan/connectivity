%{
reconstruction.OverlapGroup (computed) # Overlap densities grouping table
-> reconstruction.CellTypePair
-> reconstruction.DensityParameters
---
%}


classdef OverlapGroup < dj.Relvar & dj.AutoPopulate

	properties
		popRel = reconstruction.CellTypePair*reconstruction.DensityParameters  % !!! update the populate relation
	end

	methods(Access=protected)

		function makeTuples(self, key)
		%!!! compute missing fields for key here
			self.insert(key)
		end
	end

end