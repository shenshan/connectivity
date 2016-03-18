cell_type_names = {'L1 eNGC', 'L1 SBC', 'L23 BC', 'L23 BPC', 'L23 BTC', 'L23 ChC', 'L23 DBC', 'L23 MaC', 'L23 NGC', 'L23 Pyr','L5 BC', 'L5 DC', 'L5 HEC', 'L5 MaC', 'L5 NGC', 'L5 P', 'L5 SC'};

common_path = '/Volumes/lab/users/Shan/cell for analysis/';
for ii = 1:length(cell_type_names)
    tuple.cell_type_name = cell_type_names{ii};
    if ismember(ii, [1,2])
        tuple.layer = 'L1';
    elseif ismember(ii, 3:10)
        tuple.layer = 'L23';
    else
        tuple.layer = 'L5';
    end
    
    tuple.filepath = [common_path tuple.layer ' neurons/' tuple.cell_type_name];
    insert(reconstruction.CellType, tuple);
end