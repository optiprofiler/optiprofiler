function names = eval_fixture_select(~)
    names = {'shared', 'unloadable'};
    if strcmp(getenv('EVAL_REPORT_PLAIN_FAILURE_MODE'), 'partial'), names = {'first', 'second'}; end
end
