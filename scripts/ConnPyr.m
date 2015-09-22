
keys_from = fetch(connectivity.Cell & 'cell_layer="L5"' & 'cell_type_morph="pyr"');
keys_to = fetch(connectivity.Cell & 'cell_layer="L5"' & 'cell_type_morph="pyr"');
pairs = fetch((connectivity.ConnectMembership & 'role="from"' & keys_from)...
            * pro(connectivity.ConnectMembership & 'role="to"' & keys_to, 'cell_id->cell_id2'));
conn_pairs = fetch(connectivity.CellTestedPair & pairs & 'connected=1');

for iPair = conn_pairs'
    slice = fetch(connectivity.Slice & iPair);
    connectivity.plots.ConnSlice(slice)
end