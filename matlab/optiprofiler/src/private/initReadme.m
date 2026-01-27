function initReadme(path_readme)
%INITREADME Initialize the README file with a header.

    width = 100;
    name_width = 33;
    type_width = 8;
    gap = 2;
    try
        fid = fopen(path_readme, 'w');
        fprintf(fid, "# Content of this folder\n\n");
        % Create the header format string
        format_str = sprintf('%%-%ds%%%ds%%-%ds%%%ds%%s\\n', name_width, gap, type_width, gap);
        fprintf(fid, format_str, 'File/Folder Name', '', 'Type', '', 'Description');
        fprintf(fid, '%s\n', repmat('-', 1, width));
        fclose(fid);
    catch
    end
end
