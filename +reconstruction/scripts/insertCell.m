
cell_types = fetch(reconstruction.CellType, '*');

for iCellType = cell_types'
    
    files = dir([iCellType.filepath '/*.ASC']);
    tuple.cell_type_name = iCellType.cell_type_name;
    for ii = 1:length(files)       
        tuple.cell_id = ii;
        tuple.filename = files(ii).name;
        insert(reconstruction.Cell, tuple);
    end    
end