function addToReadme(path_readme, filename, description)
%ADDTOREADME Add an entry to the README file.

    width = 100;
    name_width = 33;
    type_width = 8;
    gap = 2;
    desc_width = width - name_width - type_width - 2 * gap;

    % Parse type and content from description
    type_str = '';
    content_str = description;
    if startsWith(description, 'File, ')
        type_str = 'File';
        content_str = description(7:end);
    elseif startsWith(description, 'Folder, ')
        type_str = 'Folder';
        content_str = description(9:end);
    end
    
    % Capitalize the first letter of the content
    if ~isempty(content_str)
        content_str(1) = upper(content_str(1));
    end

    try
        fid = fopen(path_readme, 'a');
        
        % Wrap the description using textwrap
        wrapped_desc = textwrap({content_str}, desc_width);
        
        % Manually wrap filename to ensure it breaks even without spaces
        wrapped_filename = {};
        if length(filename) <= name_width
            wrapped_filename = {filename};
        else
            current_str = filename;
            while length(current_str) > name_width
                wrapped_filename{end+1} = current_str(1:name_width);
                current_str = current_str(name_width+1:end);
            end
            if ~isempty(current_str)
                wrapped_filename{end+1} = current_str;
            end
        end
        
        n_lines = max(length(wrapped_filename), length(wrapped_desc));
        
        for i = 1:n_lines
            if i <= length(wrapped_filename)
                f_part = wrapped_filename{i};
            else
                f_part = '';
            end
            
            % Type only appears on the first line
            if i == 1
                t_part = type_str;
            else
                t_part = '';
            end
            
            if i <= length(wrapped_desc)
                d_part = wrapped_desc{i};
            else
                d_part = '';
            end
            
            format_str = sprintf('%%-%ds%%%ds%%-%ds%%%ds%%s\\n', name_width, gap, type_width, gap);
            fprintf(fid, format_str, f_part, '', t_part, '', d_part);
        end

        fclose(fid);
    catch
    end
end
