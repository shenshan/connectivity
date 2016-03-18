%{
reconstruction.CenteredTree (computed) # center the tree structure with the cell bodies
->reconstruction.Tree
-----
cellbody_coords : blob      # coordinates of the cell body, center of mass
cent_coords     : longblob  # coordinates of the cell, centered relative to the cell bodies
center_tree_ts=CURRENT_TIMESTAMP   : timestamp # automatic

%}

classdef CenteredTree < dj.Relvar & dj.AutoPopulate
	
    properties
        popRel = reconstruction.Tree
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
		
            tree = fetch(reconstruction.Tree & key, '*');
            
            cellbody = tree.node_coords(tree.node_region==2,:);
            key.cellbody_coords = mean(cellbody);
            key.cent_coords = bsxfun(@minus, tree.node_coords, key.cellbody_coords);
            
			self.insert(key)
		end
	end

end