function appendRuntimeReportNotes(results, feature, is_load, report, readme)
%APPENDRUNTIMEREPORTNOTES Keep interpretation and display changes auditable.
    fields={'fun_histories','maxcv_histories','merit_histories'};
    fid=fopen(report,'a');
    if fid<0, error('OptiProfiler:ReportOutput','Cannot append runtime notes to the report.'); end
    cleanup=onCleanup(@() fclose(fid));
    fprintf(fid,'\n## Display-only history protection\n\nFinite display values are capped at +/-1e100; raw data and scores are unchanged. Counts below describe stored history entries, including repeated tails/padding, not distinct oracle evaluations.\n');
    for k=1:numel(results)
        for j=1:numel(fields)
            values=results{k}.(fields{j});
            [clipped,nonfinite]=countDisplayEntries(values);
            fprintf(fid,'%s / %s: clipped finite entries=%d; nonfinite placeholders=%d\n', ...
                results{k}.plib,fields{j},clipped,nonfinite);
        end
    end
    if is_load
        % Feature stamps are user-customizable, not a certificate of the
        % archived oracle. Do not infer past truth semantics from their text.
        heading='Saved truth channels';
        note='Stored initial/history/output channels are retained. Replotting does not reevaluate points, repair old evaluations, or certify the original truth convention. Older quantized archives may contain inconsistent initial/output channels.';
    elseif strcmp(feature.name,'quantized')
        heading='Quantized truth';
        truth=feature.options.(FeatureOptionKey.GROUND_TRUTH.value);
        definition='original problem';
        if truth, definition='featured problem'; end
        note=sprintf('ground_truth=%d scores the %s consistently at initial/history/output points. Returned solver points are not rounded. Post-solver truth evaluation does not consume the solver budget.',truth,definition);
    else
        return;
    end
    fprintf(fid,'\n## %s\n\n%s\n',heading,note);
    addToReadme(readme,heading,note);
end

function [clipped,nonfinite]=countDisplayEntries(values)
    % Raw saved arrays may exceed 2 GB. Reporting must not allocate full-array
    % abs/logical copies merely to count entries; chunks bound temporary memory.
    clipped=0; nonfinite=0; chunk_size=1e6;
    for first=1:chunk_size:numel(values)
        chunk=values(first:min(first+chunk_size-1,numel(values)));
        finite=isfinite(chunk);
        clipped=clipped+nnz(finite & abs(chunk)>1e100);
        nonfinite=nonfinite+numel(chunk)-nnz(finite);
    end
end
