function file_name = historyPdfFileName(problem_name)
%HISTORYPDFFILENAME File name of the history plots of a problem.
% The same rule as Python's _safe_history_pdf_file_name. A report lists these
% files: a credential in a name must not reach it through the path.
    name = optiprofiler_internal.EvalReport.credentialFreeName(char(problem_name));
    file_name = [regexprep(regexprep(regexprep(strrep(name,' ','_'),'[^a-zA-Z0-9\-_]',''),'[-_]+','_'),'^[-_]+',''), '.pdf'];
end
