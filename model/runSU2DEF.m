function SU2_DEF_info=runSU2DEF...
    (mesh_filestr,cfg_filestr,dat_filestr,mesh_out_filestr,partitions,DEF_config,...
    dir_temp,run_description,remove_dump,out_logger)
% interface of SU2_DEF
% base on input filestr and DEF_config to run SU2_DEF
%
if nargin < 10
    out_logger=[];
    if nargin < 9
        remove_dump=[];
        if nargin < 8
            run_description=[];
            if nargin < 7
                dir_temp=[];
                if nargin < 6
                    DEF_config=[];
                end
            end
        end
    end
end

[str_filedir__,~]=fileparts(which('runSU2DEF.m'));

if isempty(dir_temp)
    dir_temp=fullfile(str_filedir__,'SU2_temp');
    if ~exist(dir_temp,'dir')
        mkdir(dir_temp);
    end
end

if isempty(remove_dump)
    remove_dump=false(1);
end

dir_cur=pwd();

[mesh_filename,mesh_filedir,mesh_filestr]=safeSplitFileStr(mesh_filestr);
[cfg_filename,cfg_filedir,cfg_filestr]=safeSplitFileStr(cfg_filestr);
[dat_filename,dat_filedir,dat_filestr]=safeSplitFileStr(dat_filestr);
[mesh_out_filename,mesh_out_filedir,mesh_out_filestr]=safeSplitFileStr(mesh_out_filestr);

% Config
config=Config(cfg_filestr);
config.setParameter('NUMBER_PART',partitions);

% process input
config.setParameter('DV_FILENAME',dat_filename);
config.setParameter('MESH_FILENAME',mesh_filename);
config.setParameter('MESH_OUT_FILENAME',mesh_out_filename);

% add input config
if ~isempty(DEF_config)
    field_name_list=fieldnames(DEF_config);
    for idx=1:length(field_name_list)
        field=field_name_list{idx};
        config.setParameter(field,DEF_config.(field));
    end
end

procid=feature('getpid');
if ~isempty(out_logger)
    out_logger.info(['SU2 calculating ... procid=',num2str(procid)]);
end
dir_work=fullfile(dir_temp,[run_description,'_',num2str(procid)]);

safeMakeDirs(dir_work,out_logger);
% copy mesh file to dir work
safeCopyDirs({mesh_filename},mesh_filedir,dir_work,out_logger);
% copy dat data to dir work
safeCopyDirs({dat_filename},dat_filedir,dir_work,out_logger);

% State
config.dump(fullfile(dir_work,'config_DEF.cfg'));
% state=SU2.io.State()
% out_logger.info(state)
if ~isempty(out_logger)
    out_logger.info(['mesh input file: ',mesh_filestr]);
    out_logger.info(['cfg file: ',cfg_filestr]);
    out_logger.info(['dat file: ',dat_filestr]);
    out_logger.info(['mesh input file: ',mesh_out_filestr]);
end

% run SU2 DEF
% SU2.run.DEF(config)
if ~isempty(out_logger)
    out_logger.info('begin run SU2 DEF');
end

if ispc()
    run_command='SU2_DEF config_DEF.cfg';
else
    run_command=['mpirun -np ',num2str(config.getParameter('NUMBER_PART')),' SU2_DEF config_DEF.cfg'];
end

% run command
cd(dir_work);
[status,SU2_DEF_info]=system(run_command);
cd(dir_cur);
if status ~= 0
    error_message=SU2_DEF_info;
else
    error_message=[];
    info_file=fopen(fullfile(dir_work,'SU2_DEF_info.log'),'w');
    fprintf(info_file,'%s',SU2_DEF_info);
    fclose(info_file);
    clear('info_file');
end

retry_time=0;
while (~isempty(error_message) && (retry_time < 2))
    if contains(error_message,'retrying')
        % retry again
        if ~isempty(out_logger)
            out_logger.warning('node busy,retry run SU2 DEF');
        end

        cd(dir_work);
        [status,SU2_DEF_info]=system(run_command);
        cd(dir_cur);
        if status ~= 0
            error_message=SU2_DEF_info;
        else
            error_message=[];
            info_file=fopen(fullfile(dir_work,'SU2_DEF_info.log'),'w');
            fprintf(info_file,'%s',SU2_DEF_info);
            fclose(info_file);
            clear('info_file');
        end

        retry_time=retry_time+1;
    else
        err_file=fopen(fullfile(dir_work,'SU2_DEF_err.log'),'w');
        fprintf(info_file,'%s',error_message);
        fclose(err_file);
        clear('err_file');
        if ~isempty(out_logger)
            out_logger.error([dir_work,error_message,'\n']);
        end
        error(['runSU2DEF: fatal error with SU2 DEF,proid=',num2str(procid)]);
    end
end

if ~isempty(out_logger)
    out_logger.info('end run SU2 DEF');
end

% move out
safeCopyDirs({config.getParameter('MESH_OUT_FILENAME')},dir_work,mesh_out_filedir,out_logger);

% delete temp file
if remove_dump
    if ~isempty(out_logger)
        out_logger.info('cleaning temp files...');
    end
    rmdir(dir_work,'s');
end

end
