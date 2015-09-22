function t = ConnSlice(key)
% return the connectivity table of a certain slice
%   input of the function is a key of a slice
%   output of the function is a table that shows all the connections on the
%   slice and cell type of each neuron
%   SS 2015-09-03

cells = fetch(connectivity.Cell & key,'*');
ids = [cells.cell_id];
% generate table to show connectivity
cnames = cell(1,length(cells));
rnames = cell(1,length(cells));
for ii = 1:length(cells)
    cnames{ii} = [num2str(cells(ii).cell_id) '-' cells(ii).cell_layer ' ' cells(ii).cell_type_morph];
    rnames{ii} = [num2str(cells(ii).cell_id) '-' cells(ii).cell_layer ' ' cells(ii).cell_type_morph];
end
data = zeros(length(cells),length(cells));
for ii = 1:length(cells)
    data(ii,ii) = nan;
end

conn_pairs = fetch(connectivity.CellTestedPair & key & 'connected=1');

for iPair = conn_pairs'
    cell1 = fetch(connectivity.ConnectMembership & iPair & 'role="from"');
    cell2 = fetch(connectivity.ConnectMembership & iPair & 'role="to"');
    idx1 = find(ids==cell1.cell_id);
    idx2 = find(ids==cell2.cell_id);
    data(idx1,idx2) = 1;
end

f = figure('Position',[200,200,650,250]);
t = uitable('Parent',f,'Data',data,'ColumnName',cnames,'RowName',rnames,'Position',[20,20,600,200],'ColumnEditable',false(1,length(cells)));
slice_name = fetch1(connectivity.Slice & key, 'slice_name');
title(slice_name)

