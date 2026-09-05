function results_plib = maskInvalidMerits(results_plib)
%MASKINVALIDMERITS Invalidate cached merits where saved raw evaluations are NaN.
%   Includes nested plain results. Legacy records without raw values retain
%   their merits; valid values and positive/negative infinity are unchanged.

    if isempty(results_plib)
        return;
    end
    for suffix = {'histories', 'outs', 'inits'}
        merit_key = ['merit_', suffix{1}];
        if ~isfield(results_plib, merit_key)
            continue;
        end
        merits = results_plib.(merit_key);
        if strcmp(suffix{1}, 'inits') && size(merits, 2) == 1
            for raw_key = {'fun_inits', 'maxcv_inits'}
                if isfield(results_plib, raw_key{1})
                    raw = results_plib.(raw_key{1});
                    if size(raw, 1) == size(merits, 1) && size(raw, 2) > 1
                        merits = repmat(merits, 1, size(raw, 2));
                        results_plib.(merit_key) = merits;
                        break;
                    end
                end
            end
        end
        invalid = false(size(merits));
        for prefix = {'fun_', 'maxcv_'}
            raw_key = [prefix{1}, suffix{1}];
            if isfield(results_plib, raw_key)
                invalid = invalid | isnan(results_plib.(raw_key));
            end
        end
        if any(invalid(:))
            if ~isfloat(merits)
                merits = double(merits);
            end
            merits(invalid) = Inf;
            results_plib.(merit_key) = merits;
        end
    end
    if isfield(results_plib, 'results_plib_plain') && ~isempty(results_plib.results_plib_plain)
        results_plib.results_plib_plain = maskInvalidMerits(results_plib.results_plib_plain);
    end
end
