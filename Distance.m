%{
connectivity.Distance (computed) # my newest table
-> connectivity.CellTestedPair
-----
distance:   double             # Euclidean distances between the cell pairs
dist_x  :   double             # distance in x coordinate
dist_y  :   double             # distance in y coordinate
%}

classdef Distance < dj.Relvar & dj.AutoPopulate
    properties
        popRel = pro(connectivity.CellTestedPair, connectivity.ConnectMembership*connectivity.CellPosition, 'count(*)->n') & 'n=2'
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            [x,y,role] = fetchn(connectivity.CellPosition*connectivity.ConnectMembership & key, ...
                'cell_pos_x', 'cell_pos_y', 'role');
            assert(numel(x)==2 && (all(strcmp(role,'EC')) || ~any(strcmp(role,'EC'))), 'invalid pair')
            
            % compute distances
            key.dist_x = abs(diff(x));
            key.dist_y = abs(diff(y));
            key.distance = sqrt(diff(x)^2+diff(y)^2);
            
            self.insert(key)
        end
    end
    
end