function exportDlm(from, to)
%Export connectivity and distance data to a txt file
%   each role of the output file have 
type_from = regexp(from, '(\w+) (\w+)','tokens');
type_to = regexp(to, '(\w+) (\w+)', 'tokens');

animal_adult = fetch(connectivity.Age & 'age>35');
keys_from = fetch(connectivity.Cell & animal_adult & ['cell_layer="' char(type_from{1}(1)) '"'] & ['cell_type_morph="' char(type_from{1}(2)) '"']);
keys_to = fetch(connectivity.Cell &  animal_adult & ['cell_layer="' char(type_to{1}(1)) '"'] & ['cell_type_morph="' char(type_to{1}(2)) '"']);
pairs = fetch((connectivity.ConnectMembership & 'role="from"' & keys_from)...
            * pro(connectivity.ConnectMembership & 'role="to"' & keys_to, 'cell_id->cell_id2'));
[~, connected] = fetchn(connectivity.CellTestedPair & pairs, 'pair_id', 'connected');
[dist,dist_x,dist_y] = fetchn(connectivity.Distance & pairs, 'distance', 'dist_x','dist_y');

mat = [connected, dist, dist_x, dist_y];

dlmwrite([from '_' to '.txt'],mat,'delimiter','\t','precision',3)