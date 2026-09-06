function exportPortableProfiles(curves, solver_names, options, path_stamp)
%EXPORTPORTABLEPROFILES Show exactly the computed curves without figure/JVM.
    panels = struct('title', {}, 'xlabel', {}, 'ylabel', {}, 'curves', {}, 'labels', {}, 'note', {}, 'include_zero', {});
    names = {'perf','data','log_ratio'};
    labels = {'Performance profile','Data profile','Log-ratio profile'};
    summarize = [options.(ProfileOptionKey.SUMMARIZE_PERFORMANCE_PROFILES.value), ...
        options.(ProfileOptionKey.SUMMARIZE_DATA_PROFILES.value), ...
        options.(ProfileOptionKey.SUMMARIZE_LOG_RATIO_PROFILES.value)];
    if options.(ProfileOptionKey.SEMILOGX.value)
        xlabels = {'log2(evaluation ratio)','log2(1 + evaluations/(dimension+1))','Sorted problem/run index'};
    else
        xlabels = {'Evaluation ratio','Evaluations/(dimension+1)','Sorted problem/run index'};
    end
    for i_tol = 1:numel(curves)
        for mode = {'hist','out'}
            for kind = 1:numel(names)
                record = curves{i_tol}.(mode{1});
                if ~isfield(record,names{kind}), continue; end
                source = record.(names{kind});
                if kind<3
                    source = source(:,end)'; % Same mean curves used by scores.
                    % Profiles are step functions, not linear interpolation.
                    for s=1:numel(source)
                        if isempty(source{s}), continue; end
                        x=repelem(source{s}(1,:),2); y=repelem(source{s}(2,:),2);
                        source{s}=[x(2:end);y(1:end-1)];
                    end
                end
                if all(cellfun(@isempty,source)), continue; end
                p = struct('title',sprintf('%s %s, tolerance 1e-%d',mode{1},labels{kind},i_tol), ...
                    'xlabel',xlabels{kind},'ylabel','Proportion solved', ...
                    'curves',{source},'labels',{solver_names},'note','','include_zero',kind==3);
                if kind==3, p.ylabel='log2(work solver 1 / work solver 2)'; end
                % Summary options select panels exactly as in native figures;
                % individual profiles still retain every computed curve.
                if summarize(kind) && (strcmp(mode{1},'hist') || options.(ProfileOptionKey.SUMMARIZE_OUTPUT_BASED_PROFILES.value))
                    panels(end+1)=p; %#ok<AGROW>
                end
                writeCurveSvg(fullfile(path_stamp,sprintf('%s_%s_%d.svg',names{kind},mode{1},i_tol)),p,p.title);
            end
        end
    end
    if ~isempty(panels)
        writeCurveSvg(fullfile(path_stamp,'summary.svg'),panels,'OptiProfiler profiles');
    end
    writeGraphicsIndex(path_stamp,'Native graphics or PDF export is unavailable. SVG uses the same computed curves; scores and raw data are unchanged.');
end
