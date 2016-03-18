function DistConnCorr_presaved

load('distVec.mat')

% some restrictions

idx_rel = ~isnan(connRatioVec) & nPairsVec>10 & distVec<250 & connRatioVec<1;
distVec_rel = distVec(idx_rel);
nPairsVec_rel = nPairsVec(idx_rel);
connRatioVec_rel = connRatioVec(idx_rel);



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
dist_xVec_rel = dist_xVec(idx_rel);

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
dist_yVec_rel = dist_yVec(idx_rel);
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

