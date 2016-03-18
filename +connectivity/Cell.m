%{
connectivity.Cell (manual) # my newest table
-> connectivity.Slice
cell_id                     : smallint    # 1 to 8
-----
cell_type_gene="Unknown"    : enum('Unknown','pyr','Viaat','PV','SST','VIP','5HT')  # molecular marker of the cell
cell_type_morph="Unknown"   : enum('Unknown','pyr','SBC','ENGC','MaC','BNGC','BTC','BPC','DBC','BC','ChaC','FSLP','FSLocal','ELC','DC','others')   # morphological identity of the cell
cell_type_fp="Unknown"      : enum('Unknown','FS','RS','IRS')         # spiking property of the cell
cell_layer="Unknown"        : enum('Unknown','L1','L23','L4','L5','L6') # which layer the cell body is located
cell_x                      : double                                    # position x
cell_y                      : double                                    # position y
cell_notes                  : varchar(4095)                             # other comments
cell_ts=CURRENT_TIMESTAMP   : timestamp                                 # automatic
%}

classdef Cell < dj.Relvar
end