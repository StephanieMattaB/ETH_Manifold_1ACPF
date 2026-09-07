%% GENERATE_FAIRNESS_GNE_FIGURE_PUB
%
% Publication-quality nonlinear-AC validation figure for the externally-
% supplied minimax-fair GNE x*_fair (Julia constrained-game solve,
% fairness_selected_x.csv, region CR_1): S'' baseline vs after-correction,
% all 29 physical non-slack terminals on one two-level (phase/bus) x-axis,
% security band shaded, remaining violations explicitly called out.
%
% Every layout/style knob lives in the PARAMS block below -- tune there,
% not by hunting through the plotting code. Runs unchanged in Octave or
% MATLAB.
%
% Outputs:
%   figures/fig_fairness_gne_pub.png
%   figures/fig_fairness_gne_pub.pdf
%   data/fairness_gne_pub_figure_data.mat

clear all
close all
clc

validate_fairness_gne;   % NOTE: this itself does "clear all" -- must run
                          % BEFORE the PARAMS block, not after. Produces
                          % fresh S, S', S'', x_star, corrected, Y, n,
                          % valid_mask, bus_labels, phase_labels, vmin, vmax.

%% ============================== PARAMS =================================

P = struct();

% --- canvas ---
P.fig_size_px     = [1500, 780];
P.png_dpi         = 300;
P.paper_px_per_in = 150;

% --- axes layout (figure-normalized [left bottom width height]) ---
P.main_pos      = [0.070, 0.260, 0.900, 0.560];
P.philabel_pos  = [0.070, 0.195, 0.900, 0.045];
P.buslabel_pos  = [0.070, 0.115, 0.900, 0.055];

% --- x-axis slot geometry ---
P.xpad = [0.7, 0.7];      % [left, right] padding beyond slot 1..N, data units
P.sep_color = [0.82 0.82 0.82];
P.sep_w     = 0.7;

% --- security band / reference lines ---
P.ylim               = [0.940, 1.075];
P.band_color         = [0.925 0.960 0.930];
P.band_edge_color    = [0.35 0.35 0.35];
P.band_edge_w        = 1.3;
P.refline_1p00_color = [0.70 0.70 0.70];
P.refline_1p00_w     = 0.9;

% --- series style ---
P.before_color   = [0.780 0.243 0.192];   % muted brick red -- "problem state"
P.before_marker  = '^';
P.before_msize   = 6.5;
P.before_lw      = 1.1;
P.before_ls      = ':';

P.after_color    = [0.086 0.396 0.404];   % deep teal -- "solution"
P.after_marker   = 'o';
P.after_msize    = 6.5;
P.after_lw       = 1.6;
P.after_ls       = '-';

P.resolved_edge_color = [0.086 0.396 0.404];

% --- residual-violation callouts ---
P.violation_marker_color = [0.780 0.243 0.192];
P.violation_marker_size  = 11;
P.violation_text_dy      = -0.011;   % data-units offset below the marker
P.violation_fontsize     = 8;

% --- typography ---
P.font_name           = 'Helvetica';
P.title_fontsize      = 13;
P.subtitle_fontsize   = 9.5;
P.axis_fontsize       = 9;
P.philabel_fontsize   = 8;
P.buslabel_fontsize   = 8.5;
P.legend_fontsize     = 8.5;
P.ylabel_fontsize     = 10.5;

% =========================================================================

if ~exist('figures','dir'); mkdir('figures'); end
if ~exist('data','dir');    mkdir('data');    end

phase_names = {'a','b','c'};

%% 1. Build the 29-physical-terminal slot ordering (bus-major, phase-minor)

mask_ns = valid_mask; mask_ns(1:3) = false;   % physical AND non-slack
full_idx_list = find(mask_ns);
nSlot = numel(full_idx_list);

slot_bus   = ceil(full_idx_list/3);
slot_phase = mod(full_idx_list-1,3) + 1;

v_before = Spp.ac.v(full_idx_list);
v_after  = corrected.ac.v(full_idx_list);

% bus-group boundaries, in slot-index units
bus_group_start = zeros(n,1); bus_group_end = zeros(n,1);
for b = 1:n
    idxb = find(slot_bus==b);
    if ~isempty(idxb)
        bus_group_start(b) = min(idxb);
        bus_group_end(b)   = max(idxb);
    end
end

%% 2. Identify residual violations (after-correction) for callouts

is_violated_after = (v_after < vmin) | (v_after > vmax);

%% 3. Main panel

fig = figure('visible','off', 'Position', [0 0 P.fig_size_px], 'Color','w');

ax1 = axes('Position', P.main_pos); hold(ax1,'on');
set(ax1, 'FontName', P.font_name, 'FontSize', P.axis_fontsize, 'Box','on', 'LineWidth', 0.8);

xlo = 1 - P.xpad(1); xhi = nSlot + P.xpad(2);

% security band (shaded)
patch(ax1, [xlo xhi xhi xlo], [vmin vmin vmax vmax], P.band_color, ...
    'EdgeColor','none', 'FaceAlpha', 1);
plot(ax1, [xlo xhi], [vmin vmin], '--', 'Color', P.band_edge_color, 'LineWidth', P.band_edge_w);
plot(ax1, [xlo xhi], [vmax vmax], '--', 'Color', P.band_edge_color, 'LineWidth', P.band_edge_w);
plot(ax1, [xlo xhi], [1.00 1.00], '-',  'Color', P.refline_1p00_color, 'LineWidth', P.refline_1p00_w);

% bus-group separators
for b = 1:n
    if bus_group_start(b) > 0 && bus_group_start(b) > 1
        sep = bus_group_start(b) - 0.5;
        plot(ax1, [sep sep], P.ylim, '-', 'Color', P.sep_color, 'LineWidth', P.sep_w);
    end
end

xs = 1:nSlot;
plot(ax1, xs, v_before, P.before_ls, 'Color', P.before_color, 'LineWidth', P.before_lw);
h_before = plot(ax1, xs, v_before, P.before_marker, 'MarkerSize', P.before_msize, ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', P.before_color, 'LineWidth', 1.3, 'LineStyle','none');

plot(ax1, xs, v_after, P.after_ls, 'Color', P.after_color, 'LineWidth', P.after_lw);
h_after = plot(ax1, xs, v_after, P.after_marker, 'MarkerSize', P.after_msize, ...
    'MarkerFaceColor', P.after_color, 'MarkerEdgeColor', P.after_color, 'LineStyle','none');

% residual-violation callouts
viol_idx = find(is_violated_after);
for k = viol_idx(:)'
    plot(ax1, xs(k), v_after(k), 'o', 'MarkerSize', P.violation_marker_size, ...
        'MarkerEdgeColor', P.violation_marker_color, 'LineWidth', 1.8, 'MarkerFaceColor','none');
    lbl1 = sprintf('%s.%s', bus_labels{slot_bus(k)}, phase_names{slot_phase(k)});
    lbl2 = sprintf('%.4f pu', v_after(k));
    % Two separate text() calls, not one multi-line string -- a single
    % text() with an embedded newline breaks under Octave's gnuplot
    % toolkit (emits an invalid multiplot command).
    text(ax1, xs(k), v_after(k) + P.violation_text_dy, lbl1, ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize', P.violation_fontsize, 'FontName', P.font_name, ...
        'Color', P.violation_marker_color, 'FontWeight','bold');
    text(ax1, xs(k), v_after(k) + 2*P.violation_text_dy, lbl2, ...
        'HorizontalAlignment','center', 'VerticalAlignment','top', ...
        'FontSize', P.violation_fontsize, 'FontName', P.font_name, ...
        'Color', P.violation_marker_color, 'FontWeight','bold');
end

xlim(ax1, [xlo xhi]);
ylim(ax1, P.ylim);
set(ax1, 'XTick', []);
ylabel(ax1, '|V| [p.u.]', 'FontSize', P.ylabel_fontsize, 'FontName', P.font_name);

% Built as its own dedicated strip-axes, not ax1's automatic title(),
% which under Octave's gnuplot backend renders much closer to the axes
% top than expected and collided with the subtitle/legend strip below.
ax_title = axes('Position', [0.070, 0.905, 0.900, 0.075]); %#ok<LAXES>
xlim(ax_title,[0 1]); ylim(ax_title,[0 1]); axis(ax_title,'off'); hold(ax_title,'on');
text(ax_title, 0.5, 0.3, 'Nonlinear AC Validation of the Minimax-Fair GNE x*_{fair}', ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'FontSize', P.title_fontsize, 'FontName', P.font_name, 'FontWeight','bold');

subtitle_str = sprintf(['S'''' baseline: V_{min}=%.4f pu, %d/%d violated    ->    ', ...
    'after x*_{fair}: V_{min}=%.4f pu, %d/%d violated'], ...
    Spp.rep_ac_ns.Vmin, Spp.rep_ac_ns.n_violations, Spp.rep_ac_ns.n_physical, ...
    corrected.rep_ac_ns.Vmin, corrected.rep_ac_ns.n_violations, corrected.rep_ac_ns.n_physical);
% legend() 'Position' is not respected reliably under Octave's gnuplot
% backend (it collided with in-axes text and clipped at the figure edge
% on the previous attempt). Built manually instead, in two dedicated
% blank strip-axes above the main plot, where placement is fully
% explicit and not subject to that quirk.

ax_sub = axes('Position', [0.070, 0.845, 0.500, 0.045]); %#ok<LAXES>
xlim(ax_sub,[0 1]); ylim(ax_sub,[0 1]); axis(ax_sub,'off'); hold(ax_sub,'on');
text(ax_sub, 0, 0.5, subtitle_str, 'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
    'FontSize', P.subtitle_fontsize, 'FontName', P.font_name, 'Color', [0.25 0.25 0.25]);

ax_leg = axes('Position', [0.610, 0.845, 0.360, 0.045]); %#ok<LAXES>
xlim(ax_leg,[0 1]); ylim(ax_leg,[0 1]); axis(ax_leg,'off'); hold(ax_leg,'on');
plot(ax_leg, 0.03, 0.5, P.before_marker, 'MarkerSize', P.before_msize, ...
    'MarkerFaceColor','w', 'MarkerEdgeColor', P.before_color, 'LineWidth', 1.3);
text(ax_leg, 0.09, 0.5, 'S'''' baseline', 'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
    'FontSize', P.legend_fontsize, 'FontName', P.font_name);
plot(ax_leg, 0.52, 0.5, P.after_marker, 'MarkerSize', P.after_msize, ...
    'MarkerFaceColor', P.after_color, 'MarkerEdgeColor', P.after_color);
text(ax_leg, 0.58, 0.5, 'after x*_{fair} correction', 'HorizontalAlignment','left', 'VerticalAlignment','middle', ...
    'FontSize', P.legend_fontsize, 'FontName', P.font_name);

%% 4. Two-level x-axis labels: phase letter, then bus name

ax2 = axes('Position', P.philabel_pos); xlim(ax2,[xlo xhi]); ylim(ax2,[0 1]); axis(ax2,'off'); hold(ax2,'on');
for k = 1:nSlot
    text(ax2, xs(k), 0.5, phase_names{slot_phase(k)}, 'HorizontalAlignment','center', ...
        'FontSize', P.philabel_fontsize, 'FontName', P.font_name, 'Color',[0.35 0.35 0.35]);
end

ax3 = axes('Position', P.buslabel_pos); xlim(ax3,[xlo xhi]); ylim(ax3,[0 1]); axis(ax3,'off'); hold(ax3,'on');
for b = 1:n
    if bus_group_start(b) > 0
        cx = mean([bus_group_start(b), bus_group_end(b)]);
        text(ax3, cx, 0.5, bus_labels{b}, 'HorizontalAlignment','center', ...
            'FontSize', P.buslabel_fontsize, 'FontName', P.font_name, 'FontWeight','bold');
    end
end
for b = 1:n
    if bus_group_start(b) > 0 && bus_group_start(b) > 1
        sep = bus_group_start(b) - 0.5;
        plot(ax3, [sep sep], [0 1], '-', 'Color', P.sep_color, 'LineWidth', P.sep_w);
    end
end

%% 5. Save

set(fig, 'PaperPositionMode', 'auto');
paper_w = P.fig_size_px(1)/P.paper_px_per_in;
paper_h = P.fig_size_px(2)/P.paper_px_per_in;
set(fig, 'PaperUnits', 'inches', 'PaperSize', [paper_w paper_h], 'PaperPosition', [0 0 paper_w paper_h]);
print(fig, fullfile('figures','fig_fairness_gne_pub.png'), '-dpng', sprintf('-r%d', P.png_dpi));
print(fig, fullfile('figures','fig_fairness_gne_pub.pdf'), '-dpdf');
fprintf('Saved figures/fig_fairness_gne_pub.png and .pdf\n');

save(fullfile('data','fairness_gne_pub_figure_data.mat'), ...
    'x_star','slot_bus','slot_phase','full_idx_list','v_before','v_after', ...
    'is_violated_after','bus_labels','vmin','vmax','-v7');
fprintf('Saved data/fairness_gne_pub_figure_data.mat\n');
