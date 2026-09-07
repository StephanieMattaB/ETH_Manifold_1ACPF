%% GENERATE_ANGLE_MAGNITUDE_FIGURE
%
% Reproduces the classic example_ieee13.m / example_case14.m comparison
% style -- magnitude AND angle, one-shot linear ('k*') vs nonlinear AC
% ('ko'), 3 phases x 2 columns -- for S'' specifically (the stage where the
% one-shot approximation is worst and the gap is most visible). This fills
% a real gap: every figure so far showed only |V|, never theta.
%
% Produces: figures/fig8_magnitude_angle_S2.png

clear all
close all
clc

run_ieee13_pipeline;

phnames = {'a','b','c'};

figure('visible','off','position',[0,0,900,700]);

for ph = 1:3
    rows = rw(1:n, ph);
    m = valid_mask(rows);

    v_lin = Spp.v_lin(rows); v_lin(~m) = NaN;
    v_ac  = Spp.ac.v(rows);  v_ac(~m)  = NaN;
    t_lin = Spp.t_lin(rows)/pi*180; t_lin(~m) = NaN;
    t_ac  = Spp.ac.theta(rows)/pi*180; t_ac(~m) = NaN;

    subplot(3,2,(ph-1)*2+1);
    plot(1:n, v_ac, 'ko', 1:n, v_lin, 'k*');
    xlim([0 n+1]); set(gca,'XTick',1:n,'XTickLabel',bus_labels,'FontSize',7);
    ylabel(sprintf('phase %s', phnames{ph}), 'FontWeight','bold');
    grid on;
    if ph==1; title('magnitudes |v_i| [pu]   (o = AC, * = one-shot linear)'); end

    subplot(3,2,(ph-1)*2+2);
    plot(1:n, t_ac, 'ko', 1:n, t_lin, 'k*');
    xlim([0 n+1]); set(gca,'XTick',1:n,'XTickLabel',bus_labels,'FontSize',7);
    grid on;
    if ph==1; title('angles \theta_i [deg]   (o = AC, * = one-shot linear)'); end
end

print(fullfile('figures','fig8_magnitude_angle_S2.png'), '-dpng', '-r150');
fprintf('Saved figures/fig8_magnitude_angle_S2.png\n');
