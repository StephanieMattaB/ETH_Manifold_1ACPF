function col = sequential_color(t)
%SEQUENTIAL_COLOR Perceptually-ordered single-hue colormap: pale straw
%(t=0, low deviation) -> deep red (t=1, high deviation). Deliberately not
%jet/rainbow, per the publication-figure brief (Section 7: restrained,
%single perceptual channel, low->high severity).
    t = max(0, min(1, t));
    c0 = [0.99 0.96 0.86];   % pale straw, near-white
    c1 = [0.62 0.06 0.06];   % deep red
    col = (1-t)*c0 + t*c1;
end
