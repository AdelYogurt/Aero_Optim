function [mesh_out_filestr,SU2_DEF_info]=runSU2DEF...
    (mesh_filestr,cfg_param,dat_filestr,partitions,...
    mesh_out_filestr,dir_temp,run_description,REMOVE_TEMP,out_logger)
% interface of SU2_DEF
% base on input filestr and DEF_config to run SU2_DEF
%
if nargin < 9
    out_logger=[];
    if nargin < 8
        REMOVE_TEMP=[];
        if nargin < 7
            run_description=[];
            if nargin < 6
                dir_temp=[];
                if nargin < 5
                    mesh_out_filestr=[];
                end
            end
        end
    end
end

if isempty(REMOVE_TEMP)
    REMOVE_TEMP=false(1);
end

dir_cur=pwd();

% get dir work
[str_filedir__,~]=fileparts(which('runSU2DEF.m'));
if isempty(dir_temp)
    dir_temp=fullfile(str_filedir__,'SU2_temp');
    if ~exist(dir_temp,'dir')
        mkdir(dir_temp);
    end
end
procid=feature('getpid');
dir_work=fullfile(dir_temp,[run_description,'_',num2str(procid)]);
if ~isempty(out_logger)
    out_logger.info(['SU2 calculating ... procid=',num2str(procid)]);
end

% SU2Config
if isa(cfg_param,'SU2Config')
    config=cfg_param;
    cfg_filestr=fullfile(config.config_filedir,config.config_filename);
elseif ischar(cfg_param) || isstring(cfg_param)
    cfg_filestr=cfg_param;
    [~,~,cfg_filestr]=safeSplitFileStr(cfg_filestr);
    config=SU2Config(cfg_filestr);
else
    error('runSU2CFD: error input, config is not SU2Config or cfg_filestr')
end
config.setParameter('NUMBER_PART',partitions);

% create dir work
safeMakeDirs(dir_work,out_logger);

% process mesh file
[mesh_filename,mesh_filedir,mesh_filestr]=safeSplitFileStr(mesh_filestr);
config.setParameter('MESH_FILENAME',mesh_filename);
safeCopyDirs(mesh_filename,mesh_filedir,dir_work,out_logger);

% process dat file
[dat_filename,dat_filedir,dat_filestr]=safeSplitFileStr(dat_filestr);
config.setParameter('DV_FILENAME',dat_filename);
safeCopyDirs(dat_filename,dat_filedir,dir_work,out_logger);

% process mesh out file
if isempty(mesh_out_filestr)
    mesh_out_filedir=dir_cur;
    if config.isParameter('MESH_OUT_FILENAME')
        mesh_out_filestr=fullfile(dir_cur,config.getParameter('MESH_OUT_FILENAME'));
    else
        mesh_out_filestr=fullfile(dir_cur,'mesh_out.su2');
    end
else
    [mesh_out_filename,mesh_out_filedir,mesh_out_filestr]=safeSplitFileStr(mesh_out_filestr);
    config.setParameter('MESH_OUT_FILENAME',mesh_out_filename);
end

% State
config.dump(fullfile(dir_work,'config_DEF.cfg'));
% state=SU2.io.State()
% out_logger.info(state)
if ~isempty(out_logger)
    out_logger.info(['mesh input file: ',mesh_filestr]);
    out_logger.info(['cfg file: ',cfg_filestr]);
    out_logger.info(['dat file: ',dat_filestr]);
    out_logger.info(['mesh out file: ',mesh_out_filestr]);
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
        fprintf(err_file,'%s',error_message);
        fclose(err_file);
        clear('err_file');
        if ~isempty(out_logger)
            out_logger.error([dir_work,error_message,'\n']);
        end
        error(['runSU2DEF: fatal error with SU2 DEF,proid=',num2str(procid),' error message: ',error_message]);
    end
end

if ~isempty(out_logger)
    out_logger.info('end run SU2 DEF');
end

% move mesh out
safeCopyDirs(config.getParameter('MESH_OUT_FILENAME'),dir_work,mesh_out_filedir,out_logger);

% delete temp file
if REMOVE_TEMP
    if ~isempty(out_logger)
        out_logger.info('cleaning temp files...');
    end
    rmdir(dir_work,'s');
end

end
