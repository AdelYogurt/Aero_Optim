function safeMakeDirs(dir,my_logger)
% safety make dir, error will be recorf if cannot create dir
[status,message]=mkdir(dir);
if status ~= 0 && ~isempty(my_logger)
    my_logger.error(message);
end
end
