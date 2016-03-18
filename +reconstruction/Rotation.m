%{
reconstruction.Rotation (computed) # rotate the nodes to get radially symmetric nodes
-> reconstruction.CenteredTree
-----
sample_reso   : double   # resolution of the samples, in microns
3d_samples_den  : longblob # 3d coordinates with samples on the edges for dendrites
3d_samples_axon : longblob # 3d coordinates with samples on the edges for axons
rotating_reso : double   # rotating resolution, in degs
3d_rotation_den   : longblob # 3d coordinates after rotation for dendrites
3d_rotation_axon  : longblob # 3d coordinates after rotation for dendrites
rotation_ts=CURRENT_TIMESTAMP  # automatic

%}

classdef Rotation < dj.Relvar & dj.AutoPopulate
	
    properties
        popRel = reconstruction.CenteredTree 
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            
            coords = fetch1(reconstruction.CenteredTree & key, 'cent_coords');
            [pairs, region] = fetch1(reconstruction.Tree & key, 'connected_pairs', 'node_region');
            if length(unique(node_region))<3
                return
            end
            
            sample_reso = 1;
            rotating_reso = 10;
            
            % go through all the connected node pairs
            3d_samples = 
            
			self.insert(key)
		end
	end

end