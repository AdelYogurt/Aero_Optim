function safeCopyDirs(files, old_dir, new_dir, my_logger)
% safe copy file from old dir to new dir
%
if nargin < 4
    my_logger=[];
end

if ~iscell(files)
    files={files};
end
for file_idx=1:length(files)
    file=files{file_idx};
    [status,message]=copyfile(fullfile(old_dir,file),new_dir);
    if status ~= 0
        if ~isempty(my_logger)
            my_logger.error(message);
        end
    end
end

end
