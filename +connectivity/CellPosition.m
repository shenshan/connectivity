%{
connectivity.CellPosition (computed) # get the position of the cell based on the reference axel direction (apical dendrite)
-> connectivity.Cell
-----
cell_pos_x     : double          # cell position coordinate x
cell_pos_y     : double          # cell position coordinate y
ap_direc       : double          # direction of apical dendrite, angle between with vertical -pi->pi

%}

classdef CellPosition < dj.Relvar & dj.AutoPopulate
	properties
        popRel = connectivity.Cell
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            % compute the direction of the apical dendrites
            [x1,x2,y1,y2] = fetch1(connectivity.Slice & key, 'ref1_x', 'ref2_x', 'ref1_y', 'ref2_y');
            ap_direc = atan2(x2-x1,y2-y1);
            
            % compute new coordinates of cell position
            [x,y] = fetch1(connectivity.Cell & key, 'cell_x', 'cell_y');
            [cell_pos_x, cell_pos_y] = connectivity.utils.rotateData(x,y,0,0,ap_direc,'anticlockwise');
            
            % insert key
            key.cell_pos_x = cell_pos_x;
            key.cell_pos_y = cell_pos_y;
            key.ap_direc = ap_direc;
            self.insert(key)
		end
	end

end