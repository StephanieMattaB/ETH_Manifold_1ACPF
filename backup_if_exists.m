function backup_if_exists(f)
%BACKUP_IF_EXISTS If file f already exists, rename it (timestamp suffix)
%instead of letting it be silently overwritten -- used by
%ieee13_ac_validation_final.m so re-running the validation never destroys a
%previous run's CSV output.

    if exist(f, 'file')
        [d, nm, ext] = fileparts(f);
        stamp = datestr(now, 'yyyymmdd_HHMMSS');
        movefile(f, fullfile(d, sprintf('%s_backup_%s%s', nm, stamp, ext)));
    end
end
