%% Conceptual figure

addpath('C:\Users\user\Desktop\brain_stuff\philipp\sleeppaper\src')
x = randn([500 1]);
y = colorednoise_time(0.5, 500, 1, 0, 1, 0);

subplot(2,1,1)
plot(x(1:100), 'k', 'LineWidth',3)
axis off
subplot(2,1,2)
plot(y(1:100), 'k', 'LineWidth',3)
axis off

saveas(gcf, 'C:\Users\user\Desktop\brain_stuff\te_acw_simulations\.gitignore\teplain_results\figures\conceptual_te.png')

cd C:\Users\user\Desktop\brain_stuff\te_acw_simulations
init
[v_E, time] = chaudhuri(p);
[~, ~, acf, lags] = acw(v_E(29, :), 1/p.dt);
plot(lags, acf, 'k', 'LineWidth', 3)
ax1 = gca;
yruler = ax1.YRuler;
yruler.Axle.Visible = 'off';
xruler = ax1.XRuler;
xruler.Axle.Visible = 'off'; 
xlim([0 7000])
yline(0)
yticks([0 1])
xticks([])
xlabel('Lags', 'FontSize', 15)
ylabel('ACF', 'FontSize', 15)
set(gca, 'FontSize', 15)

saveas(gcf, 'C:\Users\user\Desktop\brain_stuff\te_acw_simulations\.gitignore\teplain_results\figures\conceptual_acw.png')
%% Downsampling
cd C:\Users\user\Desktop\brain_stuff\te_acw_simulations
% init
p.w_EE = 24.3;
[v_E, time] = chaudhuri(p);
%% Downsampling figure
close all
c = 1;
irois = visual;

for i = visual([1 4 7])
    xd = rois(visual(i));
    thisroi = xd{1};
    subplot(3, 4, (1:3)+4*(c-1))
    plot(time, v_E(i, :), 'k', 'LineWidth', 2)
    xlim([52 57])
    % xticks([])
    % yticks([])
    ax1 = gca;
    yruler = ax1.YRuler;
    yruler.Visible = 'off';
    xruler = ax1.XRuler;
    xruler.Visible = 'off'; 
    ylabel(['Firing Rate (', thisroi, ')'], 'FontSize', 10)
    if c == 3, xlabel('Time', 'FontSize', 10), xruler.Label.Visible = 'on'; end
    yruler.Label.Visible = 'on';

    subplot(3, 4, 4*c)
    [acw_0, ~, acf, lags] = acw(v_E(i, :), 1/p.dt);
    plot(lags, acf, 'k', 'LineWidth', 2)
    xlim([0 350])
    ylim([-0.1 1])
    yline(0)
    % xline(acw_0)
    disp(acw_0)
    % xticks([])
    ax1 = gca;
    yruler = ax1.YRuler;
    yruler.Visible = 'off';
    xruler = ax1.XRuler;
    xruler.Visible = 'off'; 
    ylabel('ACF', 'FontSize', 10)
    if c == 3, xlabel('Lags', 'FontSize', 10), xruler.Label.Visible = 'on'; end
    yruler.Label.Visible = 'on';
    c = c + 1;
end
saveas(gcf, 'C:\Users\user\Desktop\brain_stuff\te_acw_simulations\.gitignore\teplain_results\figures\conceptual_downsampling_hier.png')
%% No hierarchy model
p.w_EE = 12.15;
[v_E, time] = chaudhuri(p);
%% Downsampling figure
close all
c = 1;
irois = visual;

for i = visual([1 4 7])
    xd = rois(visual(i));
    thisroi = xd{1};
    subplot(3, 4, (1:3)+4*(c-1))
    plot(time, v_E(i, :), 'k', 'LineWidth', 2)
    xlim([52 57])
    % xticks([])
    % yticks([])
    ax1 = gca;
    yruler = ax1.YRuler;
    yruler.Visible = 'off';
    xruler = ax1.XRuler;
    xruler.Visible = 'off'; 
    ylabel(['Firing Rate (', thisroi, ')'], 'FontSize', 10)
    if c == 3, xlabel('Time', 'FontSize', 10), xruler.Label.Visible = 'on'; end
    yruler.Label.Visible = 'on';

    subplot(3, 4, 4*c)
    [acw_0, ~, acf, lags] = acw(v_E(i, :), 1/p.dt);
    plot(lags, acf, 'k', 'LineWidth', 2)
    xlim([0 350])
    ylim([-0.1 1])
    yline(0)
    % xline(acw_0)
    disp(acw_0)
    % xticks([])
    ax1 = gca;
    yruler = ax1.YRuler;
    yruler.Visible = 'off';
    xruler = ax1.XRuler;
    xruler.Visible = 'off'; 
    ylabel('ACF', 'FontSize', 10)
    if c == 3, xlabel('Lags', 'FontSize', 10), xruler.Label.Visible = 'on'; end
    yruler.Label.Visible = 'on';
    c = c + 1;
end
saveas(gcf, 'C:\Users\user\Desktop\brain_stuff\te_acw_simulations\.gitignore\teplain_results\figures\conceptual_downsampling_nohier.png')
%% Effect of inhibition in mu_EE

cd C:\Users\user\Desktop\brain_stuff\te_acw_simulations
init
p.mu_EE = 33.7;
p.I_ext_E(1, :) = noises(1, :);
p.I_ext_E(2:29, :) = 1e-5 * randn(28, ntime);
[v_E, time, v_I] = chaudhuri(p);

subplot(2, 1, 2)
plot(time, (v_E(2, :) - mean(v_E(2, :))) / mean(v_E(2, :)), 'k')
hold on
plot(time(1:10:end), (v_I(2, 1:10:end) - mean(v_I(2, 1:10:end))) / mean(v_I(2, 1:10:end)), 'r--.')
xlim([50 51])

p.mu_EE = 32.015;
[v_E, time, v_I] = chaudhuri(p);
subplot(2, 1, 1)
plot(time, (v_E(2, :) - mean(v_E(2, :))) / mean(v_E(2, :)), 'k')
hold on
plot(time(1:10:end), (v_I(2, 1:10:end) - mean(v_I(2, 1:10:end))) / mean(v_I(2, 1:10:end)), 'r--.')
xlim([50 51])


subplot(2, 1, 2)
ylabel('Firing Rate (V2, % change)')
title('\mu_{EE} = 33.7')
legend(["Excitatory", "Inhibitory"], 'box', 'off')

subplot(2, 1, 1)
xlabel('Time (s)')
ylabel('Firing Rate (V2, % change)')
title('\mu_{EE} = 32.015')
xlabel('')
xticklabels([])


saveas(gcf, 'mu_ee_ts.png')