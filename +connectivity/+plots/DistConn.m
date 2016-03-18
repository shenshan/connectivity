function DistConn
% show the relationship between cell distance and connectivity
%   

% fetch all the pairs with their distances
pairs = fetch(connectivity.CellTestedPair);
[distances, dist_x, dist_y] = fetchn(connectivity.Distance, 'distance', 'dist_x', 'dist_y');

bins = linspace(25,500,8);
bins_x = [linspace(0,250,7),500];
bins_y = linspace(25,500,8);
conf = 0.95;
margin = 25;

z = norminv(0.5+0.5*conf);
idx = interp1(bins,1:length(bins),distances,'nearest','extrap');
idx_x = interp1(bins_x,1:length(bins_x),dist_x,'nearest','extrap');
idx_y = interp1(bins_y,1:length(bins_y),dist_y,'nearest','extrap');

connMat = zeros(1,length(bins));
errMat = zeros(1,length(bins));
num_connected_pairs = zeros(1,length(bins));
num_total_pairs = zeros(1,length(bins));
for ii = 1:length(bins)
    pairs_rel = pairs(idx==ii);
    conn =  fetchn(connectivity.CellTestedPair & pairs_rel, 'connected');
    p_conn = mean(conn);
    err = z*sqrt(p_conn*(1-p_conn)/length(conn));
    connMat(ii) = p_conn;
    errMat(ii) = err;
    num_connected_pairs(ii) = sum(conn==1);
    num_total_pairs(ii) = length(conn);
end


fig = Figure(101,'size',[60,45]); plot(bins,connMat,'ko','MarkerSize',3); hold on
errorbar(bins,connMat,errMat,'k');
ylim([0,0.3]); xlim([min(bins)-margin,max(bins)+margin]);
set(gca, 'XTick', 0:100:500)
set(gca, 'Ytick', 0:0.1:0.3)
xlabel('Intersomatic distance(\mum)');
ylabel('Connection Probability');
fig.cleanup; fig.save('plots/ConnDist.eps');

connMat_x = zeros(1,length(bins_x));
errMat_x = zeros(1,length(bins_x));
num_connected_pairs_x = zeros(1,length(bins_x));
num_total_pairs_x = zeros(1,length(bins_x));
for ii = 1:length(bins_x)
    pairs_rel = pairs(idx_x==ii);
    conn_x = fetchn(connectivity.CellTestedPair & pairs_rel, 'connected');
    p_conn_x = mean(conn_x);
    err_x = z*sqrt(p_conn_x*(1-p_conn_x)/length(conn_x));
    connMat_x(ii) = p_conn_x;
    errMat_x(ii) = err_x;
    num_connected_pairs_x(ii) = sum(conn_x==1);
    num_total_pairs_x(ii) = length(conn_x);
end

fig = Figure(102,'size',[60,45]); plot(bins_x,connMat_x,'ko', 'MarkerSize',3); hold on
errorbar(bins_x,connMat_x,errMat_x,'k');
ylim([0,0.3]); xlim([min(bins_x)-margin,max(bins_x)+margin]);
set(gca, 'XTick', 0:500:250)
set(gca, 'Ytick', 0:0.1:0.3)
xlabel('Intersomatic distance x(\mum)');
ylabel('Connection Probability');
fig.cleanup; fig.save('plots/ConnDist_x.eps');

connMat_y = zeros(1,length(bins_y));
errMat_y = zeros(1,length(bins_y));
num_connected_pairs_y = zeros(1,length(bins_y));
num_total_pairs_y = zeros(1,length(bins_y));
for ii = 1:length(bins_y)
    pairs_rel = pairs(idx_y==ii);
    conn_y = fetchn(connectivity.CellTestedPair & pairs_rel, 'connected');
    p_conn_y = mean(conn_y);
    err_y = z*sqrt(p_conn_y*(1-p_conn_y)/length(conn_y));
    connMat_y(ii) = p_conn_y;
    errMat_y(ii) = err_y;
    num_connected_pairs_y(ii) = sum(conn_y==1);
    num_total_pairs_y(ii) = length(conn_y);
end

fig = Figure(103,'size',[60,45]); plot(bins_y,connMat_y,'ko','MarkerSize',3); hold on
errorbar(bins_y,connMat_y,errMat_y, 'k');
ylim([0,0.3]); xlim([min(bins_y)-margin,max(bins_y)+margin]);
set(gca, 'XTick', 0:100:500)
set(gca, 'Ytick', 0:0.1:0.3)
xlabel('Intersomatic distance y(\mum)');
ylabel('Connection Probability');
fig.cleanup; fig.save('plots/ConnDist_y.eps');
save('All.mat', 'bins','bins_x','bins_y','num_connected_pairs','num_total_pairs','num_connected_pairs_x','num_total_pairs_x','num_connected_pairs_y','num_total_pairs_y')
