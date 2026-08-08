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
% Outputs:
%   figures/fig_voltage_topology_S_Sp_Spp.png
%   figures/fig_voltage_topology_S_Sp_Spp.pdf
%   data/voltage_topology_figure_data.mat

clear all
close all
clc

run_ieee13_pipeline;

phase_letters = {'a','b','c'};
vmin_lim = 0.95; vmax_lim = 1.05;

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

% Force consistent row-vector orientation throughout -- Octave's shape
% inference for vector-indexed-by-vector is not reliable enough to trust
% implicitly when concatenating results from several such indexings.
x_slots = x_slots(:).';
slot_full_idx = slot_full_idx(:).';
slot_phase = slot_phase(:).';

%% 2. Nonlinear-AC voltage vectors in this exact slot order

v_S   = reshape(S.ac.v(slot_full_idx),   1, []);
v_Sp  = reshape(Sp.ac.v(slot_full_idx),  1, []);
v_Spp = reshape(Spp.ac.v(slot_full_idx), 1, []);

viol_S_lo = v_S   < vmin_lim; viol_S_hi = v_S > vmax_lim;
viol_Sp_lo = v_Sp < vmin_lim; viol_Sp_hi = v_Sp > vmax_lim;
viol_Spp_lo = v_Spp < vmin_lim; viol_Spp_hi = v_Spp > vmax_lim;

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
Dmax = 10;   % shared color-scale ceiling (%), matches worst observed (~9.5%) with margin

%% 4. Topology layout (rearranged x = phase-group center; y = branch depth)

[~, edges] = ieee13_topology_layout();  % edges only (bus-index pairs); x,y recomputed below
y_bus = zeros(1,n);
y_bus([1 2 7 10]) = 0;             % trunk: 650,632,671,680
y_bus([3 4]) = 1.4;                % 633,634 (top branch)
y_bus([5 6]) = -1.4;                % 645,646 (bottom branch)
y_bus([8 9]) = 0.9;                 % 692,675
y_bus(11) = -0.9;                   % 684
y_bus(12) = -1.4;                   % 611
y_bus(13) = -0.3;                   % 652

community_buses = [5, 12, 13, 7];   % 645, 611, 652, 671

%% 5. Figure assembly

fig = figure('visible','off','position',[0,0,1700,1100]);

titles3 = {'S (baseline)', 'S'' (+IBR)', 'S'''' (+IBR, \lambda=1.8)'};
Dcols = {D(:,1), D(:,2), D(:,3)};
topo_x0 = [0.045 0.375 0.705];
topo_w  = 0.27;
topo_y0 = 0.60; topo_h = 0.36;

for si = 1:3
    ax = axes('Position', [topo_x0(si), topo_y0, topo_w, topo_h]);
    hold on;
    for e = 1:size(edges,1)
        a = edges(e,1); b = edges(e,2);
        plot([x_bus(a) x_bus(b)], [y_bus(a) y_bus(b)], '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.1);
    end
    for b = 1:n
        val = Dcols{si}(b);
        if isnan(val); col = [1 1 1]; else; col = sequential_color(val/Dmax); end
        is_comm = any(b == community_buses);
        plot(x_bus(b), y_bus(b), 'o', 'MarkerSize', 16, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.8 + 1.8*is_comm);
        text(x_bus(b), y_bus(b)-0.42, bus_labels{b}, 'HorizontalAlignment','center', ...
            'FontSize', 7, 'FontWeight','bold');
    end
    xlim([x_bus(1)-1, 30]); ylim([-2.1, 2.1]);
    axis off;
    title(titles3{si}, 'FontSize', 10);
end

% Manual sequential-color legend strip (avoids per-axes colormap issues)
axcb = axes('Position', [0.045, 0.585, 0.93, 0.02]);
ngrad = 200;
for i = 1:ngrad
    t = (i-0.5)/ngrad;
    patch([i-1 i i i-1], [0 0 1 1], sequential_color(t), 'EdgeColor', 'none');
    hold on;
end
xlim([0 ngrad]); ylim([0 1]); set(axcb,'YTick',[]);
set(axcb,'XTick', linspace(0,ngrad,6), 'XTickLabel', arrayfun(@(v) sprintf('%.0f%%',v), linspace(0,Dmax,6), 'UniformOutput', false));
xlabel(axcb, 'worst physical-phase deviation from 1.00 pu,  D_i = 100 \times max_\phi |V_{i,\phi}-1|', 'FontSize', 8);
box(axcb, 'on');

% Community-node legend marker (small inset)
axleg = axes('Position', [0.045, 0.615, 0.001, 0.001]); axis off;
annotation('ellipse', [0.045 0.955 0.012 0.018], 'FaceColor', [0.8 0.8 0.8], 'LineWidth', 2.6);
annotation('textbox', [0.062 0.945 0.30 0.03], 'String', ...
    'thick outline = future energy-community / IBR bus (645, 611, 652, 671)', ...
    'EdgeColor','none', 'FontSize', 8);

%% 6. Lower panel: phase-resolved nonlinear-AC voltage profile

ax_v = axes('Position', [0.045, 0.155, 0.93, 0.40]);
hold on;

xlo = 0.3; xhi = nPhys+0.7;
ylo = 0.895; yhi = 1.075;

% subtle security band
patch([xlo xhi xhi xlo], [vmin_lim vmin_lim vmax_lim vmax_lim], [0.94 0.97 0.94], 'EdgeColor', 'none');

% bus-group vertical separators
for b = 2:n
    if bus_group_start(b) > 1
        sep = bus_group_start(b) - 0.5;
        line([sep sep], [ylo yhi], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.7);
    end
end

% reference lines
line([xlo xhi], [1.00 1.00], 'Color', [0.6 0.6 0.6], 'LineWidth', 0.9, 'LineStyle', ':');
line([xlo xhi], [vmin_lim vmin_lim], 'Color', [0.2 0.2 0.2], 'LineWidth', 1.6);
line([xlo xhi], [vmax_lim vmax_lim], 'Color', [0.2 0.2 0.2], 'LineWidth', 1.6);

% three scenario curves -- distinct in shape, line style, AND color
hS  = plot(x_slots, v_S,  '-o', 'Color', [0.55 0.55 0.55], 'MarkerFaceColor', [0.55 0.55 0.55], ...
    'MarkerSize', 4.5, 'LineWidth', 1.1);
hSp = plot(x_slots, v_Sp, '--s', 'Color', [0.15 0.35 0.78], 'MarkerFaceColor', [0.15 0.35 0.78], ...
    'MarkerSize', 4.5, 'LineWidth', 1.4);
hSpp= plot(x_slots, v_Spp,'-^', 'Color', [0.72 0.05 0.05], 'MarkerFaceColor', [0.72 0.05 0.05], ...
    'MarkerSize', 5.5, 'LineWidth', 2.0);

% violation overlay: open black ring on any violating point, any scenario
viol_x = [x_slots(viol_S_lo | viol_S_hi), x_slots(viol_Sp_lo | viol_Sp_hi), x_slots(viol_Spp_lo | viol_Spp_hi)];
viol_y = [v_S(viol_S_lo | viol_S_hi), v_Sp(viol_Sp_lo | viol_Sp_hi), v_Spp(viol_Spp_lo | viol_Spp_hi)];
plot(viol_x, viol_y, 'o', 'MarkerSize', 11, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.3);

xlim([xlo xhi]); ylim([ylo yhi]);
set(ax_v, 'XTick', x_slots, 'XTickLabel', phase_letters(slot_phase), 'FontSize', 7.5);
ylabel('|V| [p.u.]', 'FontSize', 10);
legend([hS hSp hSpp], {'S: baseline','S'': +IBR','S'''': +IBR + congestion (\lambda=1.8)'}, ...
    'Location', 'southwest', 'FontSize', 8.5);
box(ax_v, 'on'); grid(ax_v, 'off');

% Optional min/max annotations
[vmin_val, imin] = min(v_Spp); [vmax_val, imax] = max(v_Spp);
text(x_slots(imin), vmin_val-0.010, sprintf('S'''' min: 675.c = %.4f', vmin_val), ...
    'FontSize', 7.5, 'HorizontalAlignment','center', 'Color', [0.5 0 0]);
text(x_slots(imax), vmax_val+0.010, sprintf('S'''' max: 675.b \\approx %.4f', vmax_val), ...
    'FontSize', 7.5, 'HorizontalAlignment','center', 'Color', [0.5 0 0]);

% Second label row: bus names, below the phase-letter ticks
ax_lab = axes('Position', [0.045, 0.075, 0.93, 0.055]);
xlim([xlo xhi]); ylim([0 1]); axis off; hold on;
for b = 2:n
    cx = mean([bus_group_start(b), bus_group_end(b)]);
    text(cx, 0.5, bus_labels{b}, 'HorizontalAlignment','center', 'FontSize', 8, 'FontWeight','bold');
end
for b = 2:n
    if bus_group_start(b) > 1
        sep = bus_group_start(b) - 0.5;
        line([sep sep], [0 1], 'Color', [0.85 0.85 0.85], 'LineWidth', 0.7);
    end
end

%% 7. Save

if ~exist('figures','dir'); mkdir('figures'); end
set(fig, 'PaperPositionMode', 'auto');
figpos = get(fig, 'Position');
paper_w = figpos(3)/150; paper_h = figpos(4)/150;   % inches, matching on-screen aspect
set(fig, 'PaperUnits', 'inches', 'PaperSize', [paper_w paper_h], 'PaperPosition', [0 0 paper_w paper_h]);
print(fig, fullfile('figures','fig_voltage_topology_S_Sp_Spp.png'), '-dpng', '-r300');
print(fig, fullfile('figures','fig_voltage_topology_S_Sp_Spp.pdf'), '-dpdf');
fprintf('Saved figures/fig_voltage_topology_S_Sp_Spp.png and .pdf\n');

if ~exist('data','dir'); mkdir('data'); end
save(fullfile('data','voltage_topology_figure_data.mat'), ...
    'x_slots','slot_bus','slot_phase','slot_full_idx','bus_group_start','bus_group_end', ...
    'v_S','v_Sp','v_Spp','viol_S_lo','viol_S_hi','viol_Sp_lo','viol_Sp_hi','viol_Spp_lo','viol_Spp_hi', ...
    'D','x_bus','y_bus','bus_labels','-v7');
fprintf('Saved data/voltage_topology_figure_data.mat\n');
