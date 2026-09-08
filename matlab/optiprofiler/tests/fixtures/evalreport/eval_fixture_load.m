function problem = eval_fixture_load(name)
    if strcmp(name, 'unloadable'), error('Audit:ExpectedLoadFailure', 'Deliberate public fixture load failure.'); end
    counter = getenv('EVAL_REPORT_PLAIN_FAILURE_COUNTER');
    if ~isempty(counter) && (~strcmp(getenv('EVAL_REPORT_PLAIN_FAILURE_MODE'), 'partial') || strcmp(name, 'second'))
        if isfile(counter), error('Audit:ExpectedPlainLoadFailure', 'Deliberate second-load failure.'); end
        fid = fopen(counter, 'w'); fclose(fid);
    end
    problem = Problem(struct('fun', @(x) sum(x.^2), 'x0', [1; 2], 'name', name));
end
