function write_csv(f, headers, rows)
%WRITE_CSV Minimal CSV writer (header row + cell-array data rows), used by
%ieee13_ac_validation_final.m for the three ac_validation_*.csv outputs.

    fid = fopen(f, 'w');
    fprintf(fid, '%s\n', strjoin(headers, ','));
    for r = 1:size(rows,1)
        parts = cell(1, size(rows,2));
        for c = 1:size(rows,2)
            v = rows{r,c};
            if ischar(v)
                parts{c} = v;
            elseif islogical(v)
                parts{c} = num2str(double(v));
            else
                parts{c} = sprintf('%.10g', v);
            end
        end
        fprintf(fid, '%s\n', strjoin(parts, ','));
    end
    fclose(fid);
end
