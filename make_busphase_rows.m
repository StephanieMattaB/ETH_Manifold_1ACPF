function rows = make_busphase_rows(label, v_lin, v_ac, phys_full_idx, bus_labels, phase_labels, vmin, vmax)
%MAKE_BUSPHASE_ROWS Per-terminal rows of the bus-phase AC-validation CSV
%(ieee13_ac_validation_final.m), one row per physical non-slack bus-phase
%terminal, for a single case.

    n = numel(phys_full_idx);
    rows = cell(n, 11);
    for k = 1:n
        bus = bus_labels{ceil(phys_full_idx(k)/3)};
        phase = phase_labels{mod(phys_full_idx(k)-1,3)+1};
        err = v_lin(k) - v_ac(k);
        rows(k,:) = {label, bus, phase, v_lin(k), v_ac(k), err, abs(err), ...
            v_lin(k) < vmin, v_lin(k) > vmax, v_ac(k) < vmin, v_ac(k) > vmax};
    end
end
