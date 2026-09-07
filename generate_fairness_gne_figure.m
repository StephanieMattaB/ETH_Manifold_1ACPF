%% GENERATE_FAIRNESS_GNE_FIGURE
%
% Nonlinear AC voltage-profile figure for the externally-supplied
% minimax-fair GNE x*_fair (Julia constrained-game solve,
% fairness_selected_x.csv, region CR_1), in the same visual style as
% fig5_correction_before_after.png: S'' baseline vs after-correction,
% per phase, across all 13 buses, against the [0.95,1.05] pu security
% band.
%
% Re-runs validate_fairness_gne.m fresh rather than trusting a stale
% saved .mat, so the figure always reflects the current pipeline code.

clear all
close all
clc

validate_fairness_gne;   % fresh S, S', S'', x_star, corrected (AC @ S''+x*_fair),
                          % Y, n, valid_mask, bus_labels, phase_labels, vmin, vmax

if ~exist('figures','dir'); mkdir('figures'); end

vmin_ = vmin; vmax_ = vmax;
vmask = valid_mask;

v_before = Spp.ac.v;
v_after  = corrected.ac.v;

figure('visible','off','position',[0,0,900,700]);
for ph = 1:3
    subplot(3,1,ph);
    rows = rw(1:n, ph);
    m = vmask(rows);
    xb = 1:n;
    vb = v_before(rows); vb(~m) = NaN;
    va = v_after(rows);  va(~m) = NaN;

    plot(xb, vb, 'r^-', 'MarkerFaceColor','r'); hold on;
    plot(xb, va, 'gd-', 'MarkerFaceColor','g');
    yline_(vmin_); yline_(vmax_);
    xlim([0.5, n+0.5]);
    set(gca,'XTick',1:n,'XTickLabel',bus_labels);
    ylabel(sprintf('|V| phase %s [pu]', phase_labels{ph}));
    grid on;
    if ph==1
        title('S'''' baseline vs minimax-fair-GNE-corrected x^*_{fair} (nonlinear AC)');
        legend({'S'''' baseline','after x^*_{fair}'}, 'Location', 'eastoutside');
    end
    if ph==3
        xlabel('IEEE-13 bus');
    end
end
print(fullfile('figures','fig_fairness_gne_correction.png'), '-dpng', '-r150');

fprintf('\nSaved figures/fig_fairness_gne_correction.png\n');
fprintf('  S'''' baseline : Vmin %.6f, Vmax %.6f, violations %d/%d\n', ...
    Spp.rep_ac_ns.Vmin, Spp.rep_ac_ns.Vmax, Spp.rep_ac_ns.n_violations, Spp.rep_ac_ns.n_physical);
fprintf('  after x*_fair : Vmin %.6f, Vmax %.6f, violations %d/%d\n', ...
    corrected.rep_ac_ns.Vmin, corrected.rep_ac_ns.Vmax, corrected.rep_ac_ns.n_violations, corrected.rep_ac_ns.n_physical);
