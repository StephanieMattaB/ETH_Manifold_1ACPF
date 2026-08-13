% FIG_STRESSED_VOLTAGE_PROFILE
%
% Clean phase-resolved nonlinear-AC voltage profile for the stressed
% operating point S'' only.
%
% 29 physical non-slack terminals.
% Red markers indicate violations of the 0.95--1.05 pu security band.
%
% Pure post-processing: no modification of the power-flow model.

clear all
close all
clc

run_ieee13_pipeline;

%% ============================ PARAMETERS ================================

P.fig_size_px = [1500 520];
P.png_dpi     = 300;

P.vmin = 0.95;
P.vmax = 1.05;

P.ylim = [0.895 1.075];

% Main profile
P.line_color   = [0.12 0.18 0.28];
P.line_width   = 1.5;
P.marker_size  = 5.5;

% Violations
P.viol_color   = [0.75 0.10 0.10];
P.viol_size    = 8;

% Security region
P.band_color   = [0.94 0.96 0.94];
P.limit_color  = [0.30 0.30 0.30];
P.nominal_color = [0.65 0.65 0.65];

% Bus separators
P.sep_color = [0.88 0.88 0.88];

% Typography
P.fontsize       = 9;
P.bus_fontsize   = 9;
P.ylabel_fontsize = 11;

P.show_extrema = true;

%% ====================== PHYSICAL TERMINAL ORDER =========================

phase_letters = {'a','b','c'};

x_slots = [];
slot_bus = [];
slot_phase = [];

bus_group_start = nan(1,n);
bus_group_end   = nan(1,n);

cursor = 0;

for b = 2:n       % skip slack bus 650

    first = true;

    for ph = 1:3

        fidx = rw(b,ph);

        if valid_mask(fidx)

            cursor = cursor + 1;

            x_slots(end+1)   = cursor;
            slot_bus(end+1)  = b;
            slot_phase(end+1)= ph;

            if first
                bus_group_start(b) = cursor;
                first = false;
            end

            bus_group_end(b) = cursor;

        end
    end
end

assert(numel(x_slots)==29, ...
    'Expected 29 physical non-slack terminals.');

slot_full_idx = arrayfun( ...
    @(b,p) rw(b,p), slot_bus, slot_phase);

x_slots       = x_slots(:).';
slot_bus      = slot_bus(:).';
slot_phase    = slot_phase(:).';
slot_full_idx = slot_full_idx(:).';

%% ========================== STRESSED PROFILE ============================

V = reshape(Spp.ac.v(slot_full_idx),1,[]);

viol_lo = V < P.vmin;
viol_hi = V > P.vmax;
viol    = viol_lo | viol_hi;

fprintf('\nStressed system S''''\n');
fprintf('Under-voltage terminals : %d\n',nnz(viol_lo));
fprintf('Over-voltage terminals  : %d\n',nnz(viol_hi));
fprintf('Total violations        : %d / %d\n',nnz(viol),numel(V));

%% ============================== FIGURE ==================================

fig = figure( ...
    'Color','w', ...
    'Position',[100 100 P.fig_size_px]);

ax = axes(fig);
hold(ax,'on');

xlo = 0.5;
xhi = numel(V)+0.5;

%% Security band

patch( ...
    [xlo xhi xhi xlo], ...
    [P.vmin P.vmin P.vmax P.vmax], ...
    P.band_color, ...
    'EdgeColor','none');

%% Bus separators

for b = 2:n

    if ~isnan(bus_group_start(b)) && bus_group_start(b) > 1

        xs = bus_group_start(b)-0.5;

        line( ...
            [xs xs],P.ylim, ...
            'Color',P.sep_color, ...
            'LineWidth',0.7);
    end
end

%% Reference lines

line( ...
    [xlo xhi],[1 1], ...
    'Color',P.nominal_color, ...
    'LineStyle',':', ...
    'LineWidth',1.0);

line( ...
    [xlo xhi],[P.vmin P.vmin], ...
    'Color',P.limit_color, ...
    'LineWidth',1.2);

line( ...
    [xlo xhi],[P.vmax P.vmax], ...
    'Color',P.limit_color, ...
    'LineWidth',1.2);

%% Main profile

plot( ...
    x_slots,V, ...
    '-', ...
    'Color',P.line_color, ...
    'LineWidth',P.line_width);

%% Feasible terminals

plot( ...
    x_slots(~viol),V(~viol), ...
    'o', ...
    'MarkerSize',P.marker_size, ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor',P.line_color, ...
    'LineWidth',1.1);

%% Violating terminals

plot( ...
    x_slots(viol),V(viol), ...
    'o', ...
    'MarkerSize',P.viol_size, ...
    'MarkerFaceColor',P.viol_color, ...
    'MarkerEdgeColor',P.viol_color, ...
    'LineWidth',1.1);

%% Axes

xlim([xlo xhi]);
ylim(P.ylim);

set( ...
    ax, ...
    'XTick',x_slots, ...
    'XTickLabel',phase_letters(slot_phase), ...
    'FontSize',P.fontsize, ...
    'TickDir','out', ...
    'Box','off');

ylabel('|V| [p.u.]', ...
    'FontSize',P.ylabel_fontsize);

%% ========================= EXTREMA LABELS ===============================

if P.show_extrema

    [Vmin,imin] = min(V);
    [Vmax,imax] = max(V);

    bus_min = bus_labels{slot_bus(imin)};
    ph_min  = phase_letters{slot_phase(imin)};

    bus_max = bus_labels{slot_bus(imax)};
    ph_max  = phase_letters{slot_phase(imax)};

    text( ...
        x_slots(imin),Vmin-0.008, ...
        sprintf('%s.%s  %.4f',bus_min,ph_min,Vmin), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',P.fontsize, ...
        'Color',P.viol_color);

    text( ...
        x_slots(imax),Vmax+0.008, ...
        sprintf('%s.%s  %.4f',bus_max,ph_max,Vmax), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',P.fontsize, ...
        'Color',P.viol_color);

end

%% =========================== BUS LABEL ROW ==============================

% Add bus labels just below axis using normalized text positions

yl = ylim(ax);
dy = yl(2)-yl(1);

for b = 2:n

    if isnan(bus_group_start(b))
        continue
    end

    xc = mean([bus_group_start(b),bus_group_end(b)]);

    text( ...
        xc,yl(1)-0.075*dy, ...
        bus_labels{b}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','top', ...
        'FontSize',P.bus_fontsize, ...
        'FontWeight','bold', ...
        'Clipping','off');

end

%% ================================ SAVE =================================

if ~exist('figures','dir')
    mkdir('figures');
end

set(fig,'PaperPositionMode','auto');

print( ...
    fig, ...
    fullfile('figures','fig_stressed_voltage_profile.png'), ...
    '-dpng', ...
    sprintf('-r%d',P.png_dpi));

print( ...
    fig, ...
    fullfile('figures','fig_stressed_voltage_profile.pdf'), ...
    '-dpdf');

fprintf('\nSaved stressed voltage profile.\n');