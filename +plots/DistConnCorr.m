function DistConnCorr(varargin)

% 2D plots for distances between all the cell pairs

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

distMat = zeros(length(types_from), length(types_to));
dist_xMat = zeros(length(types_from), length(types_to));
dist_yMat = zeros(length(types_from), length(types_to));
errMat = zeros(length(types_from), length(types_to));
connRatioMat = zeros(length(types_from), length(types_to));
nPairsMat = zeros(length(types_from), length(types_to));

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
        pairs = fetch(connectivity.CellTestedPair & pairs);
        [dist,dist_x,dist_y] = fetchn(connectivity.Distance & pairs,'distance','dist_x','dist_y');
        distMat(ii,jj) = mean(dist);
        dist_xMat(ii,jj) = mean(dist_x);
        dist_yMat(ii,jj) = mean(dist_y);
        errMat(ii,jj) = std(dist)/sqrt(length(dist));
        nPairsMat(ii,jj) = length(pairs);
        pairs_connected = fetch(connectivity.CellTestedPair & pairs & 'connected=1');
        connRatioMat(ii,jj) = length(pairs_connected)/length(pairs);
    end
end

distVec = distMat(:);
dist_xVec = dist_xMat(:);
dist_yVec = dist_yMat(:);
connRatioVec = connRatioMat(:);
nPairsVec = nPairsMat(:);
save('distVec.mat', 'distVec','dist_xVec','dist_yVec','connRatioVec','nPairsVec')
% some restrictions
distVec_rel = distVec(~isnan(connRatioVec));
connRatioVec_rel = connRatioVec(~isnan(connRatioVec));
nPairsVec_rel = nPairsVec(~isnan(connRatioVec));

distVec_rel = distVec_rel(nPairsVec_rel>10);
connRatioVec_rel = connRatioVec_rel(nPairsVec_rel>10);

distVec_rel = distVec_rel(distVec_rel<250);
connRatioVec_rel = connRatioVec_rel(distVec_rel<250);

distVec_rel = distVec_rel(connRatioVec_rel<0.7);
connRatioVec_rel = connRatioVec_rel(connRatioVec_rel<0.7);
% fits
p = polyfit(distVec_rel, connRatioVec_rel,1);
yfit = polyval(p,distVec_rel);
yresid = connRatioVec_rel - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(distVec_rel)-1) * var(connRatioVec_rel);
rsq = 1 - SSresid/SStotal
fig1 = Figure(101,'size',[60,45]); 
plot(distVec_rel, connRatioVec_rel,'ko','MarkerSize',3); hold on
[x,idx] = sort(distVec_rel);
plot(x, yfit(idx),'k')
set(gca, 'XTick', 0:100:400)
xlabel('Euclidean distance(\mum)')
ylabel('Connection probability')
[r,p] = corrcoef(distVec_rel, connRatioVec_rel)
fig1.cleanup; fig1.save(['plots/ConnDist_correlation_r_' num2str(r(1,2)) '.eps'])

% 
dist_xVec_rel = dist_xVec(~isnan(connRatioVec));
connRatioVec_rel = connRatioVec(~isnan(connRatioVec));
nPairsVec_rel = nPairsVec(~isnan(connRatioVec));

dist_xVec_rel = dist_xVec_rel(nPairsVec_rel>10);
connRatioVec_rel = connRatioVec_rel(nPairsVec_rel>10);

dist_xVec_rel = dist_xVec_rel(distVec_rel<250);
connRatioVec_rel = connRatioVec_rel(dist_xVec_rel<250);

dist_xVec_rel = dist_xVec_rel(connRatioVec_rel<0.7);
connRatioVec_rel = connRatioVec_rel(connRatioVec_rel<0.7);

p = polyfit(dist_xVec_rel, connRatioVec_rel,1);
yfit = polyval(p,dist_xVec_rel);

fig2 = Figure(102,'size',[60,45]); 
plot(dist_xVec_rel, connRatioVec_rel,'ko','MarkerSize',3); hold on
[x,idx] = sort(dist_xVec_rel);
plot(x, yfit(idx),'k')
xlabel('x distance(\mum)')
ylabel('Connection probability')
[r,p] = corrcoef(dist_xVec_rel, connRatioVec_rel)
fig2.cleanup; fig2.save(['plots/ConnDistX_correlation_r_' num2str(r(1,2)) '.eps'])

% 
dist_yVec_rel = dist_yVec(~isnan(connRatioVec));
connRatioVec_rel = connRatioVec(~isnan(connRatioVec));
nPairsVec_rel = nPairsVec(~isnan(connRatioVec));

dist_yVec_rel = dist_yVec_rel(nPairsVec_rel>10);
connRatioVec_rel = connRatioVec_rel(nPairsVec_rel>10);

dist_yVec_rel = dist_yVec_rel(distVec_rel<250);
connRatioVec_rel = connRatioVec_rel(dist_yVec_rel<250);

dist_yVec_rel = dist_yVec_rel(connRatioVec_rel<0.7);
connRatioVec_rel = connRatioVec_rel(connRatioVec_rel<0.7);

p = polyfit(dist_yVec_rel, connRatioVec_rel,1);
yfit = polyval(p,dist_yVec_rel);
fig3 = Figure(103,'size',[60,45]); 
plot(dist_yVec_rel, connRatioVec_rel,'ko','MarkerSize',3); hold on
[x,idx] = sort(dist_yVec_rel);
plot(x, yfit(idx),'k')
xlabel('y distance(\mum)')
ylabel('Connection probability')
[r,p] = corrcoef(dist_yVec_rel, connRatioVec_rel)
fig3.cleanup; fig3.save(['plots/ConnDistY_correlation_r_' num2str(r(1,2)) '.eps'])

