function row = make_summary_row(label, v_lin, v_ac, vmin, vmax)
%MAKE_SUMMARY_ROW One row of the case-level AC-validation summary table
%(ieee13_ac_validation_final.m), given the SAME 29 physical non-slack
%bus-phase voltages under the linear (DSO tangent) and nonlinear AC models.

    err = v_lin - v_ac;
    viol_lin = nnz(v_lin < vmin) + nnz(v_lin > vmax);
    viol_ac  = nnz(v_ac  < vmin) + nnz(v_ac  > vmax);
    row = {label, 'v_lin = DSO tangent @ S'''' AC point; v_AC = nonlinear AC replay', ...
        min(v_lin), min(v_ac), max(v_lin), max(v_ac), viol_lin, viol_ac, ...
        max(abs(err)), sqrt(mean(err.^2))};
end
