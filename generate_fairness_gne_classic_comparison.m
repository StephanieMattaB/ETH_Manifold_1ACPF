%% GENERATE_FAIRNESS_GNE_CLASSIC_COMPARISON
%
% Nonlinear-vs-linear validation figure for x*_fair, in the exact visual
% style of the original paper code's own comparison plot
% (example_ieee13.m, "%% comparison" section): one row per phase (a,b,c),
% magnitudes in the left column, angles in the right column, TRUE/reference
% values as black circles ('ko') and the model's own prediction as black
% asterisks ('k*'). No color, no shading -- deliberately the plain,
% classic academic convention used throughout the reference codebase.
%
% Unlike example_ieee13.m (which compares the published Kersting solution
% against the one-shot flat-anchored linear estimate), this compares:
%   circles (true)      = fresh, independent nonlinear AC solve at
%                          S'' + x*_fair (corrected.ac.v / .theta)
%   asterisks (predicted)= the TANGENT-PLANE linear prediction at the
%                          same point (sh.S_v/S_theta, anchored at the
%                          S'' AC operating point tp_ac) -- i.e. this is
%                          the actual linearization-accuracy check for
%                          the externally-supplied fairness GNE.
%
% Outputs:
%   figures/fig_fairness_gne_classic_comparison.png
%   figures/fig_fairness_gne_classic_comparison.pdf

clear all
close all
clc

validate_fairness_gne;   % fresh S, S', S'', x_star, corrected, sh, tp_ac,
                          % Y, n, valid_mask, bus_labels, phase_labels

if ~exist('figures','dir'); mkdir('figures'); end

%% Linear (tangent-plane) prediction at S'' + x*_fair, full non-slack vector

s_corrected_ns = [p_corrected(sh.idx.nonSlackPhase); q_corrected(sh.idx.nonSlackPhase)];
sstar_ns       = [Spp.p(sh.idx.nonSlackPhase);        Spp.q(sh.idx.nonSlackPhase)];
tstar_ns       = tp_ac.t_star(sh.idx.nonSlackPhase);
m_theta        = tstar_ns - sh.S_theta*sstar_ns;

v_pred_ns     = sh.S_v*s_corrected_ns     + sh.m_v;
theta_pred_ns = sh.S_theta*s_corrected_ns + m_theta;

% Map back onto full 3n vectors (slack rows trivially equal, since the
% slack is a fixed boundary condition, not a predicted quantity)
v_true  = corrected.ac.v;
t_true  = corrected.ac.theta;
v_pred  = v_true;   v_pred(sh.idx.nonSlackPhase)  = v_pred_ns;
t_pred  = t_true;   t_pred(sh.idx.nonSlackPhase)  = theta_pred_ns;

% Non-physical (fictitious) bus-phases: blank in both series, exactly as
% the original code's NaN gaps in v_testfeeder do
v_true(~valid_mask)  = NaN;  v_pred(~valid_mask)  = NaN;
t_true(~valid_mask)  = NaN;  t_pred(~valid_mask)  = NaN;

%% Plot: 3 rows (phase a,b,c) x 2 columns (magnitude, angle)

phnames = {'a','b','c'};

fig = figure('visible','off', 'Position', [0 0 1100 850], 'Color','w');

for ph = 1:3
    subplot(3,2,(ph-1)*2+1)
    plot(1:n, v_true(rw(1:n,ph)), 'ko', 1:n, v_pred(rw(1:n,ph)), 'k*');
    xlim([0 n+1]);
    set(gca,'XTick',[1 n],'XTickLabel',{bus_labels{1}, bus_labels{n}});
    ylabel(sprintf('phase %s', phnames{ph}), 'FontWeight','bold');
    if ph==1
        title('magnitudes  |V_i|  [pu]  --  x^*_{fair} at S''''');
    end
    if ph==3
        xlabel('bus index');
    end

    subplot(3,2,(ph-1)*2+2)
    plot(1:n, t_true(rw(1:n,ph))/pi*180, 'ko', 1:n, t_pred(rw(1:n,ph))/pi*180, 'k*');
    xlim([0 n+1]);
    set(gca,'XTick',[1 n],'XTickLabel',{bus_labels{1}, bus_labels{n}});
    if ph==1
        title('angles  \theta_i  [deg]  --  x^*_{fair} at S''''');
    end
    if ph==3
        xlabel('bus index');
    end
end

set(fig, 'PaperPositionMode', 'auto');
paper_w = 1100/150; paper_h = 850/150;
set(fig, 'PaperUnits', 'inches', 'PaperSize', [paper_w paper_h], 'PaperPosition', [0 0 paper_w paper_h]);
print(fig, fullfile('figures','fig_fairness_gne_classic_comparison.png'), '-dpng', '-r200');
print(fig, fullfile('figures','fig_fairness_gne_classic_comparison.pdf'), '-dpdf');
fprintf('Saved figures/fig_fairness_gne_classic_comparison.png and .pdf\n');

max_v_err = max(abs(v_true - v_pred));
max_t_err = max(abs(t_true - t_pred))/pi*180;
fprintf('max |v_true - v_pred| = %.3e pu\n', max_v_err);
fprintf('max |theta_true - theta_pred| = %.3e deg\n', max_t_err);
