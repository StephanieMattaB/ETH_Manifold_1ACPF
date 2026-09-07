%% GENERATE_TOPOLOGY_FIGURES
%
% Topology-map visualization of voltage severity across S, S', S'', in the
% style of a network diagram with per-node color coding (inspired by
% pandapower-style plots) -- adapted from a single-phase multi-period
% source to our three-phase, three-SCENARIO setting: each bus is colored
% by its WORST (furthest-from-1.0) physical phase voltage, nonlinear AC
% ground truth, one panel per scenario. Community nodes (645, 611, 652,
% 671) are marked with a black ring.
%
% Produces: figures/fig7_topology_map.png

clear all
close all
clc

run_ieee13_pipeline;

[xy, edges, bus_labels_topo] = ieee13_topology_layout();
community_buses = [5, 11, 13, 7];   % bus indices of 645, 611, 652, 671 (1-based, ieee13() order)

stages = {S, Sp, Spp};
titles = {'S (original)', 'S'' (+DER)', 'S'''' (congested, \lambda=1.8)'};

% Worst-phase voltage per bus, per stage.
worst_v = zeros(n, 3);
for si = 1:3
    v = stages{si}.ac.v;
    for b = 1:n
        rows = rw(b, 1:3);
        m = valid_mask(rows);
        if ~any(m); worst_v(b,si) = NaN; continue; end
        vb = v(rows); vb(~m) = NaN;
        [~, i] = max(abs(vb - 1));
        worst_v(b,si) = vb(i);
    end
end

vlo = min(worst_v(:)); vhi = max(worst_v(:));
vlo = min(vlo, 0.95); vhi = max(vhi, 1.05);

figure('visible', 'off', 'position', [0,0,1500,500]);

for si = 1:3
    subplot(1,3,si);
    hold on;

    for e = 1:size(edges,1)
        a = edges(e,1); b = edges(e,2);
        plot([xy(a,1) xy(b,1)], [xy(a,2) xy(b,2)], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.2);
    end

    for b = 1:n
        val = worst_v(b, si);
        if isnan(val)
            col = [1 1 1];
        else
            t = (val - vlo) / (vhi - vlo);
            t = max(0, min(1, t));
            col = diverging_color(t);
        end
        is_community = any(b == community_buses);
        msize = 18;
        plot(xy(b,1), xy(b,2), 'o', 'MarkerSize', msize, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1 + 2*is_community);
        text(xy(b,1), xy(b,2)-0.32, bus_labels_topo{b}, ...
            'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        if ~isnan(val)
            text(xy(b,1), xy(b,2)-0.55, sprintf('%.3f', val), ...
                'HorizontalAlignment', 'center', 'FontSize', 7);
        end
    end

    xlim([-0.8, 5.6]); ylim([-2.5, 2.4]);
    axis equal off;
    title(titles{si});
end

print(fullfile('figures','fig7_topology_map.png'), '-dpng', '-r150');
fprintf('Saved figures/fig7_topology_map.png\n');
