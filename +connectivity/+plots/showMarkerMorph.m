
function showMarkerMorph(gene)
    % show the relationship between maker and morphology
    % SS 6-19-2014
    
    % get all the cells with the maker gene
    
    cells_all = fetch(connectivity.Cell & ['cell_type_gene="' gene '"'] & 'cell_type_morph!="Unknown"');
    type_morph = unique(fetchn(connectivity.Cell & cells_all,'cell_type_morph'));
    
    perc_morph = zeros(1,length(type_morph));
    for ii = 1:length(type_morph)
        cells = fetch(connectivity.Cell & cells_all & ['cell_type_morph="' type_morph{ii} '"']);
        perc_morph(ii) = length(cells)/length(cells_all);
    end
    fig = Figure(101,'size',[90,60]); title('All'); 
    pie(perc_morph); legend(type_morph,'Location','westoutside','Orientation','vertical')
    length(cells_all)
    
    fig.cleanup; fig.save('Cell_proportion_all.eps')
    
    % L1
%     cells_all = fetch(connectivity.Cell & ['cell_type_gene="' gene '"'] & 'cell_type_morph!="Unknown"'&'cell_layer="L1"');
%     type_morph = unique(fetchn(connectivity.Cell & cells_all,'cell_type_morph'));
%     
%     perc_morph = zeros(1,length(type_morph));
%     for ii = 1:length(type_morph)
%         cells = fetch(connectivity.Cell & cells_all & ['cell_type_morph="' type_morph{ii} '"']);
%         perc_morph(ii) = length(cells)/length(cells_all);
%     end
%     fig1 = Figure(102,'size',[90,60]); title('L1');
%     pie(perc_morph); legend(type_morph,'Location','westoutside','Orientation','vertical')
%     
%     fig1.cleanup; fig1.save('Cell_proportion_L1.eps')
%     length(cells_all)
    
    % L23
    cells_all = fetch(connectivity.Cell & ['cell_type_gene="' gene '"'] & 'cell_type_morph!="Unknown"' & 'cell_layer="L23"');
    type_morph = unique(fetchn(connectivity.Cell & cells_all,'cell_type_morph'));
    
    perc_morph = zeros(1,length(type_morph));
    for ii = 1:length(type_morph)
        cells = fetch(connectivity.Cell & cells_all & ['cell_type_morph="' type_morph{ii} '"']);
        perc_morph(ii) = length(cells)/length(cells_all);
    end
    fig2 = Figure(103,'size',[90,60]); title('L23');
    pie(perc_morph); legend(type_morph,'Location','westoutside','Orientation','vertical')
    
    fig2.cleanup; fig2.save('Cell_proportion_L23');
    length(cells_all)
    
    % L5
    cells_all = fetch(connectivity.Cell & ['cell_type_gene="' gene '"'] & 'cell_type_morph!="Unknown"' & 'cell_layer="L5"');
    type_morph = unique(fetchn(connectivity.Cell & cells_all,'cell_type_morph'));
    
    perc_morph = zeros(1,length(type_morph));
    for ii = 1:length(type_morph)
        cells = fetch(connectivity.Cell & cells_all & ['cell_type_morph="' type_morph{ii} '"']);
        perc_morph(ii) = length(cells)/length(cells_all);
    end
    
    fig3 = Figure(104,'size',[90,60]); title('L5');
    pie(perc_morph); legend(type_morph,'Location','westoutside','Orientation','vertical')
    
    fig3.cleanup; fig3.save('Cell_proportion_L5');
    length(cells_all)