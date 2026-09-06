function writeGraphicsIndex(folder, reason)
%WRITEGRAPHICSINDEX Keep existing plots discoverable if aggregation fails.
    [ok,attributes]=fileattrib(folder);
    if ok, folder=attributes.Name; end
    files=[dir(fullfile(folder,'**','*.pdf'));dir(fullfile(folder,'**','*.svg'))];
    [~,order]=sort(string({files.name})); files=files(order);
    fid=fopen(fullfile(folder,'summary.html'),'w');
    if fid<0, error('OptiProfiler:OutputIndex','Cannot write output index.'); end
    cleanup=onCleanup(@() closeIfOpen(fid));
    fprintf(fid,'<!doctype html><meta charset="utf-8"><title>OptiProfiler plots</title><h1>OptiProfiler plots</h1><p>%s</p>',escapeText(reason));
    if isfile(fullfile(folder,'summary.svg'))
        fprintf(fid,'<img src="summary.svg" alt="Computed profile summary" style="max-width:100%%">');
    end
    fprintf(fid,'<ul>');
    for k=1:numel(files)
        absolute=fullfile(files(k).folder,files(k).name);
        relative=absolute(numel(folder)+2:end);
        relative=strrep(relative,filesep,'/');
        href=strrep(relative,'%','%25');
        reserved={' ','#','?','&','"','''','<','>'};
        encoded={'%20','%23','%3F','%26','%22','%27','%3C','%3E'};
        for j=1:numel(reserved), href=strrep(href,reserved{j},encoded{j}); end
        fprintf(fid,'<li><a href="%s">%s</a></li>',href,escapeText(relative));
    end
    fprintf(fid,'</ul><p>Raw numerical data and full diagnostics are retained in <a href="test_log/">test_log</a>.</p>');
    message=ferror(fid); closed=fclose(fid); clear cleanup;
    if ~isempty(message) || closed~=0
        error('OptiProfiler:OutputIndex','Could not finish the output index: %s',message);
    end
end

function text=escapeText(text)
    text=strrep(strrep(strrep(char(text),'&','&amp;'),'<','&lt;'),'>','&gt;');
    text=strrep(strrep(text,'"','&quot;'),'''','&apos;');
end

function closeIfOpen(fid)
    if ~isempty(fopen(fid)), fclose(fid); end
end
