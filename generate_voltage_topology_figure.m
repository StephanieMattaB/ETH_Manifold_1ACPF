%% GENERATE_VOLTAGE_TOPOLOGY_FIGURE
%
% Publication-quality composite figure: rearranged IEEE-13 topology
% (upper, one mini-panel per scenario, bus-level worst-phase-deviation
% color coding, shared scale) above a phase-resolved nonlinear-AC voltage
% profile (lower, all three scenarios overlaid on the same physical
% bus-phase x-axis). Uses ONLY independently solved nonlinear AC voltages
% (S.ac.v, Sp.ac.v, Spp.ac.v) -- no one-shot linear values anywhere in
% this figure. Does not touch the power-flow calculations; pure
% post-processing/plotting on already-validated pipeline results.
%
% Runs unchanged in Octave or MATLAB (no Octave-only syntax). Every
% layout/style knob lives in the PARAMS block immediately below -- tune
% there, not by hunting through the plotting code.
%
% Outputs:
%   figures/fig_voltage_topology_S_Sp_Spp.png
%   figures/fig_voltage_topology_S_Sp_Spp.pdf
%   data/voltage_topology_figure_data.mat

clear all
close all
clc

run_ieee13_pipeline;   % NOTE: this script itself does "clear all", so it must
                        % run BEFORE the PARAMS block below, not after.

%% ============================== PARAMS =================================
% Everything you're likely to want to tune lives here. Positions are
% figure-normalized [left bottom width height] in [0,1]; sizes are points
% unless noted. After changing anything, just re-run this script.

P = struct();

% --- canvas ---
P.fig_size_px   = [1700, 1100];   % [width height], on-screen/print pixels
P.png_dpi       = 300;
P.paper_px_per_in = 150;          % used to convert fig_size_px -> PaperSize (inches)

% --- upper panel: 3 topology mini-plots ---
P.topo_x0   = [0.045, 0.375, 0.705];   % left edge of each of the 3 panels
P.topo_w    = 0.27;
P.topo_y0   = 0.60;
P.topo_h    = 0.36;
P.topo_xlim_pad_left  = 1;             % data-units left of the slack node
P.topo_xlim_right     = 30;            % fixed right edge (data units)
P.topo_ylim            = [-2.1, 2.1];
P.topo_node_size       = 16;           % MarkerSize
P.topo_node_edge_w     = 0.8;          % base MarkerEdgeColor LineWidth
P.topo_community_edge_w_boost = 1.8;   % extra LineWidth for EC/IBR buses
P.topo_edge_color      = [0.65 0.65 0.65];
P.topo_edge_w          = 1.1;
P.topo_label_dy        = -0.42;        % bus-name offset below node (data units)
P.topo_label_fontsize  = 7;
P.topo_title_fontsize  = 10;
P.topo_titles = {'S (baseline)', 'S'' (+IBR)', 'S'''' (+IBR, \lambda=1.8)'};

% --- shared color scale (bus-level worst-phase-deviation metric) ---
P.Dmax = 10;   % percent; color-scale ceiling shared by all 3 topology panels
P.colorbar_pos    = [0.045, 0.585, 0.93, 0.02];   % [left bottom width height]
P.colorbar_fontsize = 8;
P.colorbar_label  = 'worst physical-phase deviation from 1.00 pu,  D_i = 100 \times max_\phi |V_{i,\phi}-1|';
P.colorbar_nticks = 6;
P.colorbar_ngrad  = 200;   % gradient resolution (patches)

% --- community/IBR legend annotation (top-left) ---
P.community_buses = [5, 12, 13, 7];   % bus indices of 645, 611, 652, 671
P.legend_ellipse_pos  = [0.045 0.955 0.012 0.018];
P.legend_textbox_pos  = [0.062 0.945 0.30 0.03];
P.legend_fontsize     = 8;
P.legend_text = 'thick outline = future energy-community / IBR bus (645, 611, 652, 671)';

% --- lower panel: voltage profile ---
P.voltage_pos     = [0.045, 0.155, 0.93, 0.40];   % [left bottom width height]
P.voltage_ylim    = [0.895, 1.075];
P.voltage_xpad    = [0.3, 0.7];        % [left_pad, right_pad] data-units beyond slot 1..29
P.vmin_lim = 0.95; P.vmax_lim = 1.05;
P.security_band_color = [0.94 0.97 0.94];
P.refline_0p95_1p05_color = [0.2 0.2 0.2];
P.refline_0p95_1p05_w     = 1.6;
P.refline_1p00_color = [0.6 0.6 0.6];
P.refline_1p00_w     = 0.9;
P.sep_color = [0.85 0.85 0.85];
P.sep_w     = 0.7;
P.tick_fontsize   = 7.5;
P.ylabel_fontsize = 10;
P.legend_series_fontsize = 8.5;
P.legend_series_location = 'southwest';

% per-scenario series style: {color, marker, linestyle, linewidth, markersize}
P.style_S   = struct('color',[0.55 0.55 0.55], 'marker','o', 'linestyle','-',  'lw',1.1, 'ms',4.5, 'label','S: baseline');
P.style_Sp  = struct('color',[0.15 0.35 0.78], 'marker','s', 'linestyle','--', 'lw',1.4, 'ms',4.5, 'label','S'': +IBR');
P.style_Spp = struct('color',[0.72 0.05 0.05], 'marker','^', 'linestyle','-',  'lw',2.0, 'ms',5.5, 'label','S'''': +IBR + congestion (\lambda=1.8)');

P.violation_ring_size = 11;
P.violation_ring_w    = 1.3;

P.show_minmax_annotations = true;   % set false if the figure looks too busy
P.minmax_fontsize = 7.5;
P.minmax_color    = [0.5 0 0];
P.minmax_dy       = 0.010;          % vertical offset of annotation text from the point

% --- second label row: bus names below the phase-letter ticks ---
P.buslabel_pos       = [0.045, 0.075, 0.93, 0.055];   % [left bottom width height]
P.buslabel_fontsize  = 8;

%% ========================================================================

phase_letters = {'a','b','c'};

%% 1. Ordered physical non-slack bus-phase x-axis (no gaps for missing phases)

x_bus = zeros(1,n);        % bus-index -> x-center (topology), all n buses incl. slack
x_slots = [];               % 1 x 29, x-coordinate of each physical non-slack slot
slot_bus = [];               % which bus (index) each slot belongs to
slot_phase = [];             % which phase (1/2/3) each slot is
bus_group_start = nan(1,n);
bus_group_end   = nan(1,n);

cursor = 0;
for b = 2:n   % skip slack (bus 1 = 650)
    first = true;
    for ph = 1:3
        fidx = rw(b, ph);
        if valid_mask(fidx)
            cursor = cursor + 1;
            x_slots(end+1) = cursor;      %#ok<AGROW>
            slot_bus(end+1) = b;          %#ok<AGROW>
            slot_phase(end+1) = ph;       %#ok<AGROW>
            if first; bus_group_start(b) = cursor; first = false; end
            bus_group_end(b) = cursor;
        end
    end
end
nPhys = numel(x_slots);   % must be 29
assert(nPhys == 29, 'Expected 29 physical non-slack terminals, got %d.', nPhys);

for b = 2:n
    x_bus(b) = mean([bus_group_start(b), bus_group_end(b)]);
end
x_bus(1) = x_bus(2) - 3;   % slack (650), placed to the left of 632

slot_full_idx = arrayfun(@(b,p) rw(b,p), slot_bus, slot_phase);

% Force consistent row-vector orientation throughout -- vector-indexed-by-
% vector shape inference is not reliable enough to trust implicitly when
% concatenating results from several such indexings (bit both Octave and,
% less often, MATLAB, so kept defensive in both).
x_slots = x_slots(:).';
slot_full_idx = slot_full_idx(:).';
slot_phase = slot_phase(:).';

%% 2. Nonlinear-AC voltage vectors in this exact slot order

v_S   = reshape(S.ac.v(slot_full_idx),   1, []);
v_Sp  = reshape(Sp.ac.v(slot_full_idx),  1, []);
v_Spp = reshape(Spp.ac.v(slot_full_idx), 1, []);

viol_S_lo = v_S   < P.vmin_lim; viol_S_hi = v_S > P.vmax_lim;
viol_Sp_lo = v_Sp < P.vmin_lim; viol_Sp_hi = v_Sp > P.vmax_lim;
viol_Spp_lo = v_Spp < P.vmin_lim; viol_Spp_hi = v_Spp > P.vmax_lim;

fprintf('S   violations: %d under / %d over\n', nnz(viol_S_lo), nnz(viol_S_hi));
fprintf('S''  violations: %d under / %d over\n', nnz(viol_Sp_lo), nnz(viol_Sp_hi));
fprintf('S'''' violations: %d under / %d over\n', nnz(viol_Spp_lo), nnz(viol_Spp_hi));

%% 3. Per-bus worst-phase-deviation metric, all 13 buses, all 3 scenarios

D = nan(n,3);   % percent, columns = [S, S', S'']
Vall = {S.ac.v, Sp.ac.v, Spp.ac.v};
for b = 1:n
    rows = rw(b,1:3);
    m = valid_mask(rows);
    if ~any(m); continue; end
    for k = 1:3
        vb = Vall{k}(rows); vb(~m) = NaN;
        D(b,k) = 100*max(abs(vb-1));
    end
end

%% 4. Topology layout (rearranged x = phase-group center; y = branch depth)

[~, edges] = ieee13_topology_layout();  % edges only (bus-index pairs); x,y recomputed above/below
y_bus = zeros(1,n);
y_bus([1 2 7 10]) = 0;             % trunk: 650,632,671,680
y_bus([3 4]) = 1.4;                % 633,634 (top branch)
y_bus([5 6]) = -1.4;                % 645,646 (bottom branch)
y_bus([8 9]) = 0.9;                 % 692,675
y_bus(11) = -0.9;                   % 684
y_bus(12) = -1.4;                   % 611
y_bus(13) = -0.3;                   % 652

%% 5. Figure assembly

fig = figure('visible','off','position',[0,0,P.fig_size_px(1),P.fig_size_px(2)]);

Dcols = {D(:,1), D(:,2), D(:,3)};

for si = 1:3
    axes('Position', [P.topo_x0(si), P.topo_y0, P.topo_w, P.topo_h]); %#ok<LAXES>
    hold on;
    for e = 1:size(edges,1)
        a = edges(e,1); b = edges(e,2);
        plot([x_bus(a) x_bus(b)], [y_bus(a) y_bus(b)], '-', 'Color', P.topo_edge_color, 'LineWidth', P.topo_edge_w);
    end
    for b = 1:n
        val = Dcols{si}(b);
        if isnan(val); col = [1 1 1]; else; col = sequential_color(val/P.Dmax); end
        is_comm = any(b == P.community_buses);
        plot(x_bus(b), y_bus(b), 'o', 'MarkerSize', P.topo_node_size, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'k', 'LineWidth', P.topo_node_edge_w + P.topo_community_edge_w_boost*is_comm);
        text(x_bus(b), y_bus(b)+P.topo_label_dy, bus_labels{b}, 'HorizontalAlignment','center', ...
            'FontSize', P.topo_label_fontsize, 'FontWeight','bold');
    end
    xlim([x_bus(1)-P.topo_xlim_pad_left, P.topo_xlim_right]); ylim(P.topo_ylim);
    axis off;
    title(P.topo_titles{si}, 'FontSize', P.topo_title_fontsize);
end

% Manual sequential-color legend strip (avoids per-axes colormap issues in
% Octave/gnuplot; in MATLAB this could instead be a real colormap+colorbar
% on a dummy image if preferred -- left as a manual strip for portability).
axcb = axes('Position', P.colorbar_pos);
for i = 1:P.colorbar_ngrad
    t = (i-0.5)/P.colorbar_ngrad;
    patch([i-1 i i i-1], [0 0 1 1], sequential_color(t), 'EdgeColor', 'none');
    hold on;
end
xlim([0 P.colorbar_ngrad]); ylim([0 1]); set(axcb,'YTick',[]);
set(axcb,'XTick', linspace(0,P.colorbar_ngrad,P.colorbar_nticks), ...
    'XTickLabel', arrayfun(@(v) sprintf('%.0f%%',v), linspace(0,P.Dmax,P.colorbar_nticks), 'UniformOutput', false));
xlabel(axcb, P.colorbar_label, 'FontSize', P.colorbar_fontsize);
box(axcb, 'on');

% Community-node legend marker (small inset)
annotation('ellipse', P.legend_ellipse_pos, 'FaceColor', [0.8 0.8 0.8], 'LineWidth', 2.6);
annotation('textbox', P.legend_textbox_pos, 'String', P.legend_text, ...
    'EdgeColor','none', 'FontSize', P.legend_fontsize);

%% 6. Lower panel: phase-resolved nonlinear-AC voltage profile

ax_v = axes('Position', P.voltage_pos);
hold on;

xlo = 1 - P.voltage_xpad(1); xhi = nPhys + P.voltage_xpad(2);
ylo = P.voltage_ylim(1); yhi = P.voltage_ylim(2);

% subtle security band
patch([xlo xhi xhi xlo], [P.vmin_lim P.vmin_lim P.vmax_lim P.vmax_lim], P.security_band_color, 'EdgeColor', 'none');

% bus-group vertical separators
for b = 2:n
    if bus_group_start(b) > 1
        sep = bus_group_start(b) - 0.5;
        line([sep sep], [ylo yhi], 'Color', P.sep_color, 'LineWidth', P.sep_w);
    end
end

% reference lines
line([xlo xhi], [1.00 1.00], 'Color', P.refline_1p00_color, 'LineWidth', P.refline_1p00_w, 'LineStyle', ':');
line([xlo xhi], [P.vmin_lim P.vmin_lim], 'Color', P.refline_0p95_1p05_color, 'LineWidth', P.refline_0p95_1p05_w);
line([xlo xhi], [P.vmax_lim P.vmax_lim], 'Color', P.refline_0p95_1p05_color, 'LineWidth', P.refline_0p95_1p05_w);

% three scenario curves -- distinct in shape, line style, AND color
hS  = plot(x_slots, v_S,  ['-' P.style_S.marker],  'Color', P.style_S.color,  'LineStyle', P.style_S.linestyle, ...
    'MarkerFaceColor', P.style_S.color,  'MarkerSize', P.style_S.ms,  'LineWidth', P.style_S.lw);
hSp = plot(x_slots, v_Sp, ['-' P.style_Sp.marker], 'Color', P.style_Sp.color, 'LineStyle', P.style_Sp.linestyle, ...
    'MarkerFaceColor', P.style_Sp.color, 'MarkerSize', P.style_Sp.ms, 'LineWidth', P.style_Sp.lw);
hSpp= plot(x_slots, v_Spp,['-' P.style_Spp.marker],'Color', P.style_Spp.color,'LineStyle', P.style_Spp.linestyle, ...
    'MarkerFaceColor', P.style_Spp.color,'MarkerSize', P.style_Spp.ms,'LineWidth', P.style_Spp.lw);

% violation overlay: open black ring on any violating point, any scenario
viol_x = [x_slots(viol_S_lo | viol_S_hi), x_slots(viol_Sp_lo | viol_Sp_hi), x_slots(viol_Spp_lo | viol_Spp_hi)];
viol_y = [v_S(viol_S_lo | viol_S_hi), v_Sp(viol_Sp_lo | viol_Sp_hi), v_Spp(viol_Spp_lo | viol_Spp_hi)];
plot(viol_x, viol_y, 'o', 'MarkerSize', P.violation_ring_size, 'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'none', 'LineWidth', P.violation_ring_w);

xlim([xlo xhi]); ylim([ylo yhi]);
set(ax_v, 'XTick', x_slots, 'XTickLabel', phase_letters(slot_phase), 'FontSize', P.tick_fontsize);
ylabel('|V| [p.u.]', 'FontSize', P.ylabel_fontsize);
legend([hS hSp hSpp], {P.style_S.label, P.style_Sp.label, P.style_Spp.label}, ...
    'Location', P.legend_series_location, 'FontSize', P.legend_series_fontsize);
box(ax_v, 'on'); grid(ax_v, 'off');

if P.show_minmax_annotations
    [vmin_val, imin] = min(v_Spp); [vmax_val, imax] = max(v_Spp);
    minb = bus_labels{slot_bus(imin)}; minph = phase_letters{slot_phase(imin)};
    maxb = bus_labels{slot_bus(imax)}; maxph = phase_letters{slot_phase(imax)};
    text(x_slots(imin), vmin_val-P.minmax_dy, sprintf('S'''' min: %s.%s = %.4f', minb, minph, vmin_val), ...
        'FontSize', P.minmax_fontsize, 'HorizontalAlignment','center', 'Color', P.minmax_color);
    text(x_slots(imax), vmax_val+P.minmax_dy, sprintf('S'''' max: %s.%s \\approx %.4f', maxb, maxph, vmax_val), ...
        'FontSize', P.minmax_fontsize, 'HorizontalAlignment','center', 'Color', P.minmax_color);
end

% Second label row: bus names, below the phase-letter ticks
axes('Position', P.buslabel_pos); %#ok<LAXES>
xlim([xlo xhi]); ylim([0 1]); axis off; hold on;
for b = 2:n
    cx = mean([bus_group_start(b), bus_group_end(b)]);
    text(cx, 0.5, bus_labels{b}, 'HorizontalAlignment','center', 'FontSize', P.buslabel_fontsize, 'FontWeight','bold');
end
for b = 2:n
    if bus_group_start(b) > 1
        sep = bus_group_start(b) - 0.5;
        line([sep sep], [0 1], 'Color', P.sep_color, 'LineWidth', P.sep_w);
    end
end

%% 7. Save

if ~exist('figures','dir'); mkdir('figures'); end
set(fig, 'PaperPositionMode', 'auto');
paper_w = P.fig_size_px(1)/P.paper_px_per_in;
paper_h = P.fig_size_px(2)/P.paper_px_per_in;
set(fig, 'PaperUnits', 'inches', 'PaperSize', [paper_w paper_h], 'PaperPosition', [0 0 paper_w paper_h]);
print(fig, fullfile('figures','fig_voltage_topology_S_Sp_Spp.png'), '-dpng', sprintf('-r%d', P.png_dpi));
print(fig, fullfile('figures','fig_voltage_topology_S_Sp_Spp.pdf'), '-dpdf');
fprintf('Saved figures/fig_voltage_topology_S_Sp_Spp.png and .pdf\n');

if ~exist('data','dir'); mkdir('data'); end
save(fullfile('data','voltage_topology_figure_data.mat'), ...
    'x_slots','slot_bus','slot_phase','slot_full_idx','bus_group_start','bus_group_end', ...
    'v_S','v_Sp','v_Spp','viol_S_lo','viol_S_hi','viol_Sp_lo','viol_Sp_hi','viol_Spp_lo','viol_Spp_hi', ...
    'D','x_bus','y_bus','bus_labels','-v7');
fprintf('Saved data/voltage_topology_figure_data.mat\n');
