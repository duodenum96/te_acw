%% Eta analysis
clear, clc, close all
cd /BICNAS2/ycatal/te_acw/modeling/.gitignore/teplain_results
netw = "eta";
tes = load(netw + "netw_task_te_plain.mat");
te  = tes.te;
acws_rest = load("acws_" + netw + "netw_rest.mat");
acw_rest = acws_rest.acw0;

% ACW: dims: rois x sims x manips
% tes: dims: rois x rois x sims x manips

nroi = size(acw_rest, 1);
nsim = size(acw_rest, 2);
nmanip = size(acw_rest, 3);

nd = zeros(nroi, nsim, nmanip);
for m = 1:nmanip
    for i = 1:nsim
        for j = 1:nroi
            nd(j, i, m) = sum(te(j, :, i, m)) + sum(te(:, j, i, m));
        end
    end
end
manips = [0.8 0.82 0.84 0.86 0.88 0.9 0.92 0.94 0.96 0.98 1];
nmanip = length(manips);
teslope = zeros(1, nmanip);
acwslope = zeros(1, nmanip);
rois = repmat((1:nroi)', nsim, 1);
for i = 1:nmanip
    slope = polyfit(rois, squeeze(acw_rest(:, :, i)), 1);
    acwslope(i) = slope(1);

    slope = polyfit(rois, squeeze(nd(:, :, i)), 1);
    teslope(i) = slope(1);
end

[rho, p] = corr(acwslope', teslope', 'type','Spearman');

%% 
cmap = cbrewer('qual', 'Accent', nmanip); cmap = cmap ./ max(cmap);

colormap(cmap)
scatter(acwslope, teslope, 70, flipud(cmap), 'filled', 'MarkerEdgeColor', 'k')
cb = colorbar('Ticks', 1/22:1/11:1, 'TickLabels', 0.68 * flip(manips));
ylabel(cb, "\eta", "Rotation", 0, "FontSize", 14)
l = lsline;
l.Color = 'r';
xlabel('Slope of linear fit to ROIs vs ACW')
ylabel('Slope of linear fit to ROIs vs total TE')
annotation('textbox',[0.642785714285711,0.803333333333335,0.14614285714286,0.121904761904772], ...
    'String',"\rho = " + num2str(rho, 3) + " p < 0.001",'EdgeColor','none')
grid on
saveas(gcf, 'figures\eta_additional.png')