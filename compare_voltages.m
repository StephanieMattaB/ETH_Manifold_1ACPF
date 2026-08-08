function cmp = compare_voltages(v_lin, v_ac, valid_mask, bus_labels, phase_labels)
%COMPARE_VOLTAGES Manifold/linear vs nonlinear-AC voltage-magnitude error,
%restricted to physically connected bus-phases.
%
%   cmp = COMPARE_VOLTAGES(v_lin, v_ac, valid_mask, bus_labels, phase_labels)

    idx_phys = find(valid_mask);
    err = v_lin(idx_phys) - v_ac(idx_phys);

    [cmp.max_error, i_worst] = max(abs(err));
    cmp.rms_error = sqrt(mean(err.^2));
    cmp.error_phys = err;
    cmp.idx_phys = idx_phys;

    [bus, ph] = deal_terminal(idx_phys(i_worst), bus_labels, phase_labels);
    cmp.worst_bus = bus;
    cmp.worst_phase = ph;
    cmp.worst_v_lin = v_lin(idx_phys(i_worst));
    cmp.worst_v_ac  = v_ac(idx_phys(i_worst));
end


function [bus_name, phase_name] = deal_terminal(terminal_idx, bus_labels, phase_labels)
    bus_idx   = ceil(terminal_idx/3);
    phase_idx = mod(terminal_idx - 1, 3) + 1;
    bus_name   = bus_labels{bus_idx};
    phase_name = phase_labels{phase_idx};
end
