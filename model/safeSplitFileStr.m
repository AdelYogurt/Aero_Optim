function [filename, filedir, filestr]=safeSplitFileStr(filestr)
% safety split file str to file name and file dir
%
[filedir,name,ext]=fileparts(filestr);
filename=[name,ext];
if ~exist(filedir,"dir")
    filedir=fullfile(pwd(),filedir);
    filestr=fullfile(pwd(),filestr);
end
end