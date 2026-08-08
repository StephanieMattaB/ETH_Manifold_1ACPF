function col = diverging_color(t)
%DIVERGING_COLOR blue (low, t=0) -> green (mid, t=0.5) -> red (high, t=1).
%Simple 3-stop RGB interpolation, t in [0,1].
    t = max(0, min(1, t));
    if t < 0.5
        f = t/0.5;
        col = (1-f)*[0.10 0.30 0.85] + f*[0.15 0.70 0.20];
    else
        f = (t-0.5)/0.5;
        col = (1-f)*[0.15 0.70 0.20] + f*[0.80 0.10 0.10];
    end
end
