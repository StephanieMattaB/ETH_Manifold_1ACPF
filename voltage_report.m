function rep = voltage_report(v, valid_mask, bus_labels, phase_labels, vmin, vmax)
%VOLTAGE_REPORT Common voltage-security diagnostics for a 3n x 1 voltage
%magnitude vector, restricted to physically connected bus-phases.
%
%   rep = VOLTAGE_REPORT(v, valid_mask, bus_labels, phase_labels, vmin, vmax)
%
%   v            : 3n x 1 voltage magnitudes [p.u.].
%   valid_mask   : 3n x 1 logical, true at physically connected terminals.
%   bus_labels   : 1 x n cell array of bus names, in ieee13() bus order.
%   phase_labels : 1 x 3 cell array, e.g. {'a','b','c'}.
%   vmin, vmax   : voltage-security bounds (e.g. 0.95, 1.05).
%
%   Records exactly the diagnostics requested at every stage of the DSO
%   pipeline: Vmin/Vmax and their bus-phase, the full physical-terminal
%   voltage list, under/over-voltage counts, and the identity of every
%   violated bus-phase.

    idx_phys = find(valid_mask);
    v_phys   = v(idx_phys);

    [rep.Vmin, i_min] = min(v_phys);
    [rep.Vmax, i_max] = max(v_phys);

    [rep.Vmin_bus, rep.Vmin_phase] = terminal_name(idx_phys(i_min), bus_labels, phase_labels);
    [rep.Vmax_bus, rep.Vmax_phase] = terminal_name(idx_phys(i_max), bus_labels, phase_labels);

    under = v_phys < vmin;
    over  = v_phys > vmax;

    rep.n_undervoltage = nnz(under);
    rep.n_overvoltage  = nnz(over);
    rep.n_violations   = rep.n_undervoltage + rep.n_overvoltage;
    rep.n_physical     = numel(idx_phys);

    rep.v_phys     = v_phys;
    rep.idx_phys   = idx_phys;
    rep.vmin       = vmin;
    rep.vmax       = vmax;

    viol_idx = idx_phys([find(under); find(over)]);
    rep.violated_terminals = cell(numel(viol_idx), 1);
    for k = 1:numel(viol_idx)
        [b, ph] = terminal_name(viol_idx(k), bus_labels, phase_labels);
        rep.violated_terminals{k} = sprintf('%s.%s (%.4f pu)', b, ph, v(viol_idx(k)));
    end
end


function [bus_name, phase_name] = terminal_name(terminal_idx, bus_labels, phase_labels)
    bus_idx   = ceil(terminal_idx/3);
    phase_idx = mod(terminal_idx - 1, 3) + 1;
    bus_name   = bus_labels{bus_idx};
    phase_name = phase_labels{phase_idx};
end
