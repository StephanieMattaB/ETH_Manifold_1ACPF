function s = ternary(cond, a, b)
%TERNARY Small helper: s = a if cond else b.
    if cond
        s = a;
    else
        s = b;
    end
end
