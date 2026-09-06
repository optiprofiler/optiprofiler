function exportPortableHistory(file, mode, problem_name, problem_dim, solver_names, histories, inits, n_eval, labels, options)
%EXPORTPORTABLEHISTORY The native mean/errorband curves without graphics.
    panels=struct('title',{},'xlabel',{},'ylabel',{},'curves',{},'labels',{},'note',{});
    rows = false;
    if strcmp(mode,'cummin'), rows=true; elseif strcmp(mode,'combined'), rows=[false,true]; end
    for cumulative=rows
        for channel=1:numel(histories)
            [values,note]=processHistYaxes(histories{channel},inits{channel});
            shift=computeHistoryYShift(values,options);
            [x,means,lower,upper,n_runs]=prepareHistoryPlotData(values,cumulative,shift,n_eval,options);
            is_log=false;
            for s=1:numel(means)
                is_log=is_log || (any(means{s}) && any(diff(means{s})));
            end
            lines={}; names={};
            for solver=1:numel(means)
                series={means{solver}}; suffix={''};
                if n_runs>1
                    series={means{solver},lower{solver},upper{solver}};
                    suffix={' mean',' lower band',' upper band'};
                end
                for j=1:numel(series)
                    y=series{j};
                    if is_log
                        y(y<=0)=NaN; y=log10(y);
                    end
                    lines{end+1}=[x{solver}/(problem_dim+1);y]; %#ok<AGROW>
                    names{end+1}=[solver_names{solver},suffix{j}]; %#ok<AGROW>
                end
            end
            title=sprintf('%s: %s',problem_name,labels{channel});
            if cumulative, title=[title,' (cumulative minimum)']; end
            ylabel=sprintf('Value + shift %.4g',shift);
            if is_log, ylabel=sprintf('log10(value + shift %.4g)',shift); end
            panels(end+1)=struct('title',title,'xlabel','Evaluations/(dimension+1)', ...
                'ylabel',ylabel,'curves',{lines},'labels',{names},'note',note); %#ok<AGROW>
        end
    end
    writeCurveSvg(file,panels,problem_name);
end
