function [mesh_out_filestr,dir_work]=runSU2DEF...
    (cfg_param,partitions,dir_temp,run_desc,REMOVE_TEMP,out_logger)
% interface of SU2_DEF
% base on input cfg_param and partitions to run SU2_DEF
%
% input:
% cfg_param: Cfg_filestr or class of 'SU2Config'.
% partitions: Processes number of parallel run SU2_DEF
%
% output:
% mesh_out_filestr: Filestr of output mesh
% SU2_DEF_info: Shell output of SU2_DEF.
%
if nargin < 6
    out_logger=[];
    if nargin < 5
        REMOVE_TEMP=[];
        if nargin < 4
            run_desc=[];
            if nargin < 3
                dir_temp=[];
            end
        end
    end
end

if isempty(REMOVE_TEMP)
    REMOVE_TEMP=false(1);
end

dir_cur=pwd();

% SU2Config
if isa(cfg_param,'SU2Config')
    config=cfg_param;
    cfg_filestr=fullfile(config.config_filedir,config.config_filename);
elseif ischar(cfg_param) || isstring(cfg_param)
    cfg_filestr=cfg_param;
    [~,~,cfg_filestr]=safeSplitFileStr(cfg_filestr);
    config=SU2Config(cfg_filestr);
else
    error('runSU2DEF: error input, config is not SU2Config or cfg_filestr')
end
config.setParameter('NUMBER_PART',partitions);

% get dir work
if isempty(dir_temp)
    dir_temp=fullfile(pwd(),'SU2_temp');
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
if ~isempty(out_logger)
    out_logger.info(['SU2 calculating ... procid=',num2str(procid)]);
end

% create dir work
safeMakeDirs(dir_work,out_logger);

% process mesh file
if ~config.isParameter('MESH_FILENAME')
    error('runSU2DEF: config lack MESH_FILENAME define');
end
mesh_filestr=config.getParameter('MESH_FILENAME');
[mesh_filename,mesh_filedir,mesh_filestr]=safeSplitFileStr(mesh_filestr);
if ~exist(mesh_filestr,'file')
    error('runSU2DEF: mesh file do not exist')
end
config.setParameter('MESH_FILENAME',mesh_filename);
if contains(mesh_filename,'cgns','IgnoreCase',true)
    config.setParameter('MESH_FORMAT','CGNS');
elseif contains(mesh_filename,'su2','IgnoreCase',true)
    config.setParameter('MESH_FORMAT','SU2');
else
    error('runSU2DEF: unsupport mesh type');
end
safeCopyDirs(mesh_filename,mesh_filedir,dir_work,out_logger);

% process dat file
if config.isParameter('DV_FILENAME')
    dat_filestr=config.getParameter('DV_FILENAME');
    [dat_filename,dat_filedir,~]=safeSplitFileStr(dat_filestr);
    config.setParameter('DV_FILENAME',dat_filename);
    safeCopyDirs(dat_filename,dat_filedir,dir_work,out_logger);
end

% process mesh out file
if config.isParameter('MESH_OUT_FILENAME')
    mesh_out_filestr=fullfile(dir_cur,config.getParameter('MESH_OUT_FILENAME'));
else
    mesh_out_filestr=fullfile(dir_cur,'mesh_out.su2');
end
[mesh_out_filename,mesh_out_filedir,mesh_out_filestr]=safeSplitFileStr(mesh_out_filestr);
config.setParameter('MESH_OUT_FILENAME',mesh_out_filename);

% State
config.dump(fullfile(dir_work,'config_DEF.cfg'));
% state=SU2.io.State()
% out_logger.info(state)
if ~isempty(out_logger)
    out_logger.info(['mesh input file: ',mesh_filestr]);
    out_logger.info(['cfg file: ',cfg_filestr]);
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
    run_command=['mpirun -np ',num2str(partitions),' SU2_DEF config_DEF.cfg'];
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
            out_logger.error([dir_work,error_message]);
        end
        error(['runSU2DEF: fatal error with SU2 DEF,proid=',num2str(procid)]);
    end
end

if ~isempty(out_logger)
    out_logger.info('end run SU2 DEF');
end

% move mesh out
safeCopyDirs(mesh_out_filename,dir_work,mesh_out_filedir,out_logger);

% delete temp file
if REMOVE_TEMP
    if ~isempty(out_logger)
        out_logger.info(['cleaning temp directory and files: ',dir_work]);
    end
    rmdir(dir_work,'s');
end

end
