function ConnSpec(varargin)

% 2D plots for connectivity between all the cell pairs

if (length(varargin) == 1 && strcmp(lowercase(varargin{1}),'all')) || isempty(varargin)
    % get all the cell types in the database
    types_from = unique(fetchn(connectivity.Cell & ...
    'cell_layer not in ("L4","L6","Unknown")' & 'cell_type_morph!="Unknown"',...
    'CONCAT(cell_layer," ",cell_type_morph)->cell_type'));
    types_to = types_from;
else
    idx1 = find(strcmp(varargin,'from')==1);
    idx2 = find(strcmp(varargin, 'to')==1);
    types_from = varargin{idx1+1:idx2-1};
    types_to = varargin{idx2+1:end};
end

connNumMat = zeros(length(types_from), length(types_to));
allNumMat = zeros(length(types_from), length(types_to));
connRatioMat = zeros(length(types_from), length(types_to));

for ii = 1:length(types_from)
    for jj = 1:length(types_to)
        type_from = types_from{ii};
        type_to = types_to{jj};
        type_from = regexp(type_from, '(\w+) (\w+)','tokens');
        type_to = regexp(type_to, '(\w+) (\w+)', 'tokens');
        keys_from = fetch(connectivity.Cell & ['cell_layer="' char(type_from{1}(1)) '"'] & ['cell_type_morph="' char(type_from{1}(2)) '"']);
        keys_to = fetch(connectivity.Cell & ['cell_layer="' char(type_to{1}(1)) '"'] & ['cell_type_morph="' char(type_to{1}(2)) '"']);
        pairs = fetch((connectivity.ConnectMembership & 'role="from"' & keys_from)...
            * pro(connectivity.ConnectMembership & 'role="to"' & keys_to, 'cell_id->cell_id2'));
        pairs_all = fetch(connectivity.CellTestedPair & pairs);
        pairs_connected = fetch(connectivity.CellTestedPair & pairs & 'connected=1');
        connNumMat(ii,jj) = length(pairs_connected);
        allNumMat(ii,jj) = length(pairs_all);
        connRatioMat(ii,jj) = length(pairs_connected)/length(pairs_all);
    end
end

figure; set(gcf,'Position',[50,50,1000,800]);
imagesc(1:length(types_from), 1:length(types_to),connRatioMat');
set(gca,'XTick',1:length(types_from),'XTickLabel', types_from, 'YTick',1:length(types_to), 'YTickLabel', types_to); colormap(jet)

for ii = 1:length(types_from)
    for jj = 1:length(types_to)
        text(ii-0.3,jj,[num2str(connNumMat(ii,jj)) '/' num2str(allNumMat(ii,jj))])
    end
end