function [fluent_out,fluent_info]=runFluentCFD...
    (mesh_filestr,jou_filestr,partitions,fluent_dir,dimension,out_filename,...
    dir_temp,run_desc)
% interface of fluent CFD
%
if nargin < 8
    run_desc=[];
    if nargin < 7
        dir_temp=[];
    end
end

dir_cur=pwd();

% get dir work
[str_filedir__,~]=fileparts(which('runFluentCFD.m'));
if isempty(dir_temp)
    dir_temp=fullfile(str_filedir__,'Fluent_temp');
    if ~exist(dir_temp,'dir')
        mkdir(dir_temp);
    end
end
procid=feature('getpid');
if ~isempty(run_desc)
    dir_work=fullfile(dir_temp,run_desc);
else
    dir_work=fullfile(dir_temp,['PID_',num2str(procid)]);
end

% create dir work
safeMakeDirs(dir_work);

% process mesh file
[mesh_filename,mesh_filedir,~]=safeSplitFileStr(mesh_filestr);
safeCopyDirs(mesh_filename,mesh_filedir,dir_work);

% process jou file
[jou_filename,jou_filedir,~]=safeSplitFileStr(jou_filestr);
safeCopyDirs(jou_filename,jou_filedir,dir_work);

% build command
if ispc()
    fluent_filestr=fullfile(fluent_dir,'fluent.exe');
else
    fluent_filestr=fullfile(fluent_dir,'fluent');
end
run_command=['"',fluent_filestr,'" ',num2str(dimension,'%d'),'ddp ',...
    '-wait ','-g ',num2str(partitions,'-t%d'),' ','-i ',jou_filename];

% run command
cd(dir_work);

% check if exist old file
delete *.h5
delete *.trn
delete *.out

[status,fluent_info]=system(run_command);
cd(dir_cur);

if status ~= 0
    error_message=fluent_info;

    err_file=fopen(fullfile(dir_work,'fluent_err.log'),'w');
    fprintf(err_file,'%s',error_message);
    fclose(err_file);
    clear('err_file');
    error(['runFluentCFD: fatal error with fluent,proid=',num2str(procid),' error message: ',error_message]);
else
    info_file=fopen(fullfile(dir_work,'fluent_info.log'),'w');
    fprintf(info_file,'%s',fluent_info);
    fclose(info_file);
    clear('info_file');
end

% obtain result
out_filestr=fullfile(dir_work,out_filename);
fluent_out=importdata(out_filestr);

end