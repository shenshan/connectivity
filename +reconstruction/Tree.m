%{
reconstruction.Tree (imported) # read neurolucida file as a tree structure
-> reconstruction.Cell
-----
node_num            : int       # number of nodes of this tree
node_coords         : longblob  # x,y,z coordinates of each node, nNodes x 3
connected_pairs     : longblob  # connected pairs of nodes, nConnectedPairs x 2
node_region         : longblob  # the region a node belongs to, 1 is Axon, 2 is CellBody, 3 is Dendrite
tree_ts=CURRENT_TIMESTAMP   : timestamp # automatic
%}

classdef Tree < dj.Relvar & dj.AutoPopulate
	
    properties
        popRel = reconstruction.Cell
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            tic
            filepath = fetch1(reconstruction.CellType & key, 'filepath');
            filename = fetch1(reconstruction.Cell & key, 'filename');
            cbuf = [filepath '/' filename];
            tree = neurolucida_tree(cbuf);
            key.node_num = length(tree.X);
            key.node_coords = [tree.X,tree.Y,tree.Z];
            
            if length(tree.rnames)~=3 || sum(strcmp(tree.rnames, {'Axon' 'CellBody' 'Dendrite'}))~=3
            	
                if ismember('Axon', tree.rnames)
                    idx = strcmp(tree.rnames, 'Axon');
                    id = find(idx==1);
                    tree.R(tree.R==id) = 1;
                end
                
                if ismember('CellBody', tree.rnames)
                    idx = strcmp(tree.rnames, 'CellBody');
                    id = find(idx==1);
                    tree.R(tree.R==id) = 2;
                end
                
                if ismember('Dendrite', tree.rnames)
                    idx = strcmp(tree.rnames, 'Dendrite');
                    id = find(idx==1);
                    tree.R(tree.R==id) = 3;
                end
            end
            key.node_region = tree.R;
            [i,j] = find(tree.dA);
            key.connected_pairs = [i,j];
            toc
			self.insert(key)
		end
	end

end