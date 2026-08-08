function q_cap = ieee13_capacitors(n)
%IEEE13_CAPACITORS Fixed shunt-capacitor reactive injections of the original
%IEEE-13 benchmark, as a 3n x 1 vector in the same [bus-major, a/b/c] layout
%used throughout this repository (see rw.m).
%
%   q_cap = IEEE13_CAPACITORS(n)
%
%   These are exactly the two capacitor banks that ieee13.m folds into its
%   net-consumption vector q by subtraction (s3(9,:) -= 200i on all three
%   phases of bus 9 = "675"; s3(12,3) -= 100i on phase c of bus 12 = "611").
%   q_cap holds the POSITIVE var values that were subtracted, i.e. q_cap(k)
%   is the capacitor's injection magnitude at bus-phase k (0 where none).
%
%   Reconstructing them here, once, avoids re-deriving this same hard-coded
%   pair of numbers independently in every script that needs to separate
%   "load" from "net consumption" (congestion scaling must scale load only
%   and leave capacitor banks fixed, per the DSO-pipeline stress convention).
%
%   Units: per-unit on Sbase = 5e6 VA (matches ieee13.m).

    Sbase = 5e6;

    q_cap = zeros(3*n, 1);
    q_cap(rw(9))    = [200; 200; 200] * 1e3 / Sbase;   % bus 675, phases a,b,c
    q_cap(rw(12,3)) = 100e3 / Sbase;                    % bus 611, phase c
end
