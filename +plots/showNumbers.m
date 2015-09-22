function showNumbers

types = unique(fetchn(connectivity.Cell & ...
    'cell_layer not in ("L4","L6","Unknown")' & 'cell_type_morph!="Unknown"',...
    'CONCAT(cell_layer," ",cell_type_morph)->cell_type'));

% cell number
cell_numbers_vec = cell(2,length(types));

for ii = 1:length(types)
    cell_numbers_vec{1,ii} = types{ii};
    type = regexp(types{ii}, '(\w+) (\w+)','tokens');
    cell_numbers_vec{2,ii} = length(fetch(connectivity.Cell & ['cell_layer="' char(type{1}(1)) '"'] & ['cell_type_morph="' char(type{1}(2)) '"']));
end

xlwrite('CellNumber.xls', cell_numbers_vec');


% pairs number
pairs_num_mat = cell(length(types)+1,length(types)+1);
conn_pairs_num_mat = cell(size(pairs_num_mat));

for ii = 1:size(pairs_num_mat,1)
    for jj = 1:size(pairs_num_mat,2)
        if ii == 1 && jj == 1
            pairs_num_mat{ii,jj} = '';
            conn_pairs_num_mat{ii,jj} = '';
        elseif ii == 1 && jj ~= 1
            pairs_num_mat{1,jj} = types{jj-1};
            conn_pairs_num_mat{ii,jj} = types{jj-1};
        elseif ii ~= 1 && jj == 1
            pairs_num_mat{ii,1} = types{ii-1};
            conn_pairs_num_mat{ii,jj} = types{ii-1};
        else
            type_from = types{ii-1};
            type_to = types{jj-1};
            type_from = regexp(type_from, '(\w+) (\w+)','tokens');
            type_to = regexp(type_to, '(\w+) (\w+)', 'tokens');
            keys_from = fetch(connectivity.Cell & ['cell_layer="' char(type_from{1}(1)) '"'] & ['cell_type_morph="' char(type_from{1}(2)) '"']);
            keys_to = fetch(connectivity.Cell & ['cell_layer="' char(type_to{1}(1)) '"'] & ['cell_type_morph="' char(type_to{1}(2)) '"']);
            pairs = fetch((connectivity.ConnectMembership & 'role="from"' & keys_from)...
                * pro(connectivity.ConnectMembership & 'role="to"' & keys_to, 'cell_id->cell_id2'));
            pairs = fetch(connectivity.CellTestedPair & pairs);
            conn_pairs = fetch(connectivity.CellTestedPair & pairs & 'connected=1');
            pairs_num_mat{ii,jj} = length(pairs);
            conn_pairs_num_mat{ii,jj} = length(conn_pairs);
        end
            
    end
end

xlwrite('PairsNumber.xls',pairs_num_mat,1);
xlwrite('PairsNumber.xls',conn_pairs_num_mat,2);
