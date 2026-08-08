%% GENERATE_REPORT_FIGURES
%
% Produces the figure set for the DSO-pipeline research write-up:
%   fig1_voltage_profile.png       - per-phase V profile, S vs S' vs S'' (AC)
%   fig2_violations.png            - violation counts per stage (bar)
%   fig3_linearization_error.png   - one-shot linear vs AC error per stage
%   fig4_sensitivity_heatmap.png   - |S_v| community sensitivity heatmap
%   fig5_correction_before_after.png - S'' baseline vs min-effort-corrected AC
%   (fig6 reuses the already-generated data/io_map_accuracy.png)
%
% Loads results from the already-run, already-saved pipeline stages rather
% than re-deriving anything -- pure post-processing/plotting.

clear all
close all
clc

run_ieee13_pipeline;   % fresh S, S', S'', sh, Y, n, valid_mask, bus_labels, ...
fr = load(fullfile('data','feasibility_results.mat'));

if ~exist('figures','dir'); mkdir('figures'); end

phase_labels = {'a','b','c'};
vmin_ = 0.95; vmax_ = 1.05;

%% fig1: voltage profile per phase, three stages, nonlinear AC

v_S = S.ac.v; v_Sp = Sp.ac.v; v_Spp = Spp.ac.v;
vmask = valid_mask;

figure('visible','off','position',[0,0,900,700]);
for ph = 1:3
    subplot(3,1,ph);
    rows = rw(1:n, ph);
    m = vmask(rows);
    xb = 1:n;
    vS_ph = v_S(rows); vS_ph(~m) = NaN;
    vSp_ph = v_Sp(rows); vSp_ph(~m) = NaN;
    vSpp_ph = v_Spp(rows); vSpp_ph(~m) = NaN;

    plot(xb, vS_ph, 'ko-', 'MarkerFaceColor','k'); hold on;
    plot(xb, vSp_ph, 'bs-', 'MarkerFaceColor','b');
    plot(xb, vSpp_ph, 'r^-', 'MarkerFaceColor','r');
    yline_(vmin_); yline_(vmax_);
    xlim([0.5, n+0.5]);
    set(gca,'XTick',1:n,'XTickLabel',bus_labels);
    ylabel(sprintf('|V| phase %s [pu]', phase_labels{ph}));
    grid on;
    if ph==1
        title('Nonlinear AC voltage profile: S vs S'' vs S'''' (lambda=1.8)');
        legend('S','S''','S''''','Location','southwest');
    end
    if ph==3
        xlabel('IEEE-13 bus');
    end
end
print(fullfile('figures','fig1_voltage_profile.png'), '-dpng', '-r150');

%% fig2: violation counts per stage (non-slack physical terminals)

under = [S.rep_ac_ns.n_undervoltage, Sp.rep_ac_ns.n_undervoltage, Spp.rep_ac_ns.n_undervoltage];
over  = [S.rep_ac_ns.n_overvoltage,  Sp.rep_ac_ns.n_overvoltage,  Spp.rep_ac_ns.n_overvoltage];

figure('visible','off','position',[0,0,500,400]);
bar([under; over]', 'stacked');
set(gca,'XTickLabel',{'S','S''','S'''''});
ylabel('non-slack violation count (of 29)');
legend('undervoltage','overvoltage','Location','northwest');
title('Voltage-security violations per stage (nonlinear AC)');
grid on;
print(fullfile('figures','fig2_violations.png'), '-dpng', '-r150');

%% fig3: one-shot linearization error vs nonlinear AC, per stage

maxerr = [S.cmp.max_error, Sp.cmp.max_error, Spp.cmp.max_error];
rmserr = [S.cmp.rms_error, Sp.cmp.rms_error, Spp.cmp.rms_error];

figure('visible','off','position',[0,0,600,420]);
bar([maxerr; rmserr]');
set(gca,'XTickLabel',{'S','S''','S'''''}, 'YScale','log');
ylim([1e-4, 1]);
ylabel('voltage error vs nonlinear AC [pu]');
legend({'max error','RMS error'}, 'Location', 'eastoutside');
title('One-shot manifold linearization error (grows under stress)');
grid on;
print(fullfile('figures','fig3_linearization_error.png'), '-dpng', '-r150');

%% fig4: |S_v| community sensitivity heatmap (physical rows x 14 community cols)

Sv_F = sh.A_F_sh(1:sh.nPhys, :);   % top half of A_F_sh IS S_v in community column order

physicalFull = sh.idx.nonSlackPhase(sh.valid_ns);
row_labels = cell(sh.nPhys,1);
for i = 1:sh.nPhys
    bus_idx = ceil(physicalFull(i)/3);
    ph_idx  = mod(physicalFull(i)-1,3)+1;
    row_labels{i} = sprintf('%s.%s', bus_labels{bus_idx}, phase_labels{ph_idx});
end
col_labels = {'p645b','p645c','p611c','p652a','p671a','p671b','p671c', ...
              'q645b','q645c','q611c','q652a','q671a','q671b','q671c'};

figure('visible','off','position',[0,0,700,900]);
imagesc(abs(Sv_F));
colorbar;
set(gca,'XTick',1:14,'XTickLabel',col_labels,'XTickLabelRotation',90);
set(gca,'YTick',1:sh.nPhys,'YTickLabel',row_labels,'FontSize',7);
title('|S_v| sensitivity: community [p,q] -> physical terminal voltage [pu/pu]');
print(fullfile('figures','fig4_sensitivity_heatmap.png'), '-dpng', '-r150');

%% fig5: S'' baseline vs minimum-effort-corrected AC voltage profile

v_before = Spp.ac.v; v_after = fr.corrected.ac.v;

figure('visible','off','position',[0,0,900,700]);
for ph = 1:3
    subplot(3,1,ph);
    rows = rw(1:n, ph);
    m = vmask(rows);
    xb = 1:n;
    vb = v_before(rows); vb(~m) = NaN;
    va = v_after(rows); va(~m) = NaN;

    plot(xb, vb, 'r^-', 'MarkerFaceColor','r'); hold on;
    plot(xb, va, 'gd-', 'MarkerFaceColor','g');
    yline_(vmin_); yline_(vmax_);
    xlim([0.5, n+0.5]);
    set(gca,'XTick',1:n,'XTickLabel',bus_labels);
    ylabel(sprintf('|V| phase %s [pu]', phase_labels{ph}));
    grid on;
    if ph==1
        title('S'''' baseline vs minimum-effort-corrected (nonlinear AC)');
        legend({'S'''' baseline','after correction'}, 'Location', 'eastoutside');
    end
    if ph==3
        xlabel('IEEE-13 bus');
    end
end
print(fullfile('figures','fig5_correction_before_after.png'), '-dpng', '-r150');

fprintf('Saved 5 figures to figures/\n');
