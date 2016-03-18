
key.cell_type_name = 'L5 MaC';
key.cell_id = 1;

coords = fetch1(reconstruction.CenteredTree & key, 'cent_coords');
[pairs, region] = fetch1(reconstruction.Tree & key, 'connected_pairs', 'node_region');
 
% figure; hold on
% 
% scatter3(coords(region==1,1),coords(region==1,2),coords(region==1,3),'r.')
% scatter3(coords(region==2,1),coords(region==2,2),coords(region==2,3),'k.')
% scatter3(coords(region==3,1),coords(region==3,2),coords(region==3,3),'g.')
% 
% xlabel('x'); ylabel('y'); zlabel('z')

% separate dendrites and axons

x = pairs(:,1); y = pairs(:,2);
idx_den = find(region==3);
idx_axon = find(region==1);
pairs_den = pairs(ismember(x,idx_den)&ismember(y,idx_den),:);
pairs_axon = pairs(ismember(x,idx_axon)&ismember(y,idx_axon),:);

% make samples on edges of dendrites and axons
