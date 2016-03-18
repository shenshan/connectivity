function DistHist(type_from, type_to)
%DISTHIST plots the histogram of the distance of certain cell types
%   type_from and type_to specify the cell types pre and post synaptic
%   SS 15-09-06
from = type_from; to = type_to;
type_from = regexp(type_from, '(\w+) (\w+)','tokens');
type_to = regexp(type_to, '(\w+) (\w+)', 'tokens');
keys_from = fetch(connectivity.Cell & ['cell_layer="' char(type_from{1}(1)) '"'] & ['cell_type_morph="' char(type_from{1}(2)) '"']);
keys_to = fetch(connectivity.Cell & ['cell_layer="' char(type_to{1}(1)) '"'] & ['cell_type_morph="' char(type_to{1}(2)) '"']);
pairs = fetch((connectivity.ConnectMembership & 'role="from"' & keys_from)...
    * pro(connectivity.ConnectMembership & 'role="to"' & keys_to, 'cell_id->cell_id2'));
pairs = fetch(connectivity.CellTestedPair & pairs);
[dist,dist_x,dist_y] = fetchn(connectivity.Distance & pairs,'distance','dist_x','dist_y');

fig = Figure(101,'size',[80,50]);

[n,x] = hist(dist(dist<250));
dist_rel = dist(dist<250);
mean(dist_rel)
std(dist_rel)/sqrt(length(dist_rel))
h = bar(x,n);
set(h,'facecolor','w','barwidth',1)
xlabel('Euclidean distance (\mum)')
ylabel('Counts')

fig.cleanup
fig.save([from '_' to '_Euclidean_hist.eps'])