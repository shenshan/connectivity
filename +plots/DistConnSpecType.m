function DistConnSpecType(type_from, type_to, layer_from, layer_to, tosave)
% show the relationship between cell distance and connectivity, with
% specific layer conditions
%   
if ~exist('tosave','var')
    tosave=0;
end
% fetch all the pairs with their distances
layer_from_keys = fetch(connectivity.Cell & ['cell_layer="' layer_from '"'] & ['cell_type_morph="' type_from '"']);
layer_to_keys = fetch(connectivity.Cell & ['cell_layer="' layer_to '"'] & ['cell_type_morph="' type_to '"']);
pairs = fetch((connectivity.ConnectMembership & 'role="from"' & layer_from_keys)...
    * pro(connectivity.ConnectMembership & 'role="to"' & layer_to_keys, 'cell_id->cell_id2'));
pairs = fetch(connectivity.CellTestedPair & pairs);
[distances, dist_x, dist_y] = fetchn(connectivity.Distance & pairs, 'distance', 'dist_x', 'dist_y');

bins = 50:20:200;
bins_x = 50:20:200;
bins_y = 50:20:200;
conf = 0.95;

z = norminv(0.5+0.5*conf);
idx = interp1(bins,1:length(bins),distances,'nearest','extrap');
idx_x = interp1(bins_x,1:length(bins_x),dist_x,'nearest','extrap');
idx_y = interp1(bins_y,1:length(bins_y),dist_y,'nearest','extrap');

connMat = zeros(1,length(bins));
errMat = zeros(1,length(bins));
for ii = 1:length(bins)
    pairs_rel = pairs(idx==ii);
    conn =  fetchn(connectivity.CellTestedPair & pairs_rel, 'connected');
    p_conn = mean(conn);
    err = z*sqrt(p_conn*(1-p_conn)/length(conn));
    connMat(ii) = p_conn;
    errMat(ii) = err;
end

fig = Figure(101,'size',[80,60]); plot(bins,connMat,'o'); hold on
errorbar(bins,connMat,errMat);
ylim([0,0.3]); xlim([min(bins)-50,max(bins)+50]);
xlabel('Intersomatic distance(um)');
ylabel('Connection Probability');
fig.cleanup; 

if tosave
    fig.save(['ConnDist_' layer_from '_' layer_to '.eps']);
end

connMat_x = zeros(1,length(bins_x));
errMat_x = zeros(1,length(bins_x));
for ii = 1:length(bins_x)
    pairs_rel = pairs(idx_x==ii);
    conn_x = fetchn(connectivity.CellTestedPair & pairs_rel, 'connected');
    p_conn_x = mean(conn_x);
    err_x = z*sqrt(p_conn_x*(1-p_conn_x)/length(conn_x));
    connMat_x(ii) = p_conn_x;
    errMat_x(ii) = err_x;
end

fig = Figure(102,'size',[80,60]); plot(bins_x,connMat_x,'o'); hold on
errorbar(bins_x,connMat_x,errMat_x);
ylim([0,0.3]); xlim([min(bins_x)-50,max(bins_x)+50]);
xlabel('Intersomatic distance x(um)');
ylabel('Connection Probability');
fig.cleanup;

if tosave
    fig.save(['ConnDist_x_' layer_from '_' layer_to '.eps']);
end

connMat_y = zeros(1,length(bins_y));
errMat_y = zeros(1,length(bins_y));
for ii = 1:length(bins_y)
    pairs_rel = pairs(idx_y==ii);
    conn_y = fetchn(connectivity.CellTestedPair & pairs_rel, 'connected');
    p_conn_y = mean(conn_y);
    err_y = z*sqrt(p_conn_y*(1-p_conn_y)/length(conn_y));
    connMat_y(ii) = p_conn_y;
    errMat_y(ii) = err_y;
end

fig = Figure(103,'size',[80,60]); plot(bins_y,connMat_y,'o'); hold on
errorbar(bins_y,connMat_y,errMat_y);
ylim([0,0.3]); xlim([min(bins_y)-50,max(bins_y)+50]);
xlabel('Intersomatic distance y(um)');
ylabel('Connection Probability');
fig.cleanup;

if tosave
    fig.save(['ConnDist_y_' layer_from '_' layer_to '.eps']);
end