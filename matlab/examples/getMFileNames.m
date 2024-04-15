function fileList = getMFileNames(folderPath, excludeList)

    files = dir(fullfile(folderPath, '*.m'));
    fileList = {};

    for i = 1:length(files)
        filename = files(i).name;
        filename = filename(1:end-2);

        if ~ismember(filename, excludeList)
            fileList{end+1} = filename;
        end
    end
end