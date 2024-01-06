function [SU2_data,SU2_history,SU2_surface,SU2_CFD_info]=runSU2CFD...
    (cfg_param,partitions,dir_temp,run_desc,REMOVE_TEMP,out_logger)
% interface of SU2_CFD
% base on input cfg_param and partitions to run SU2_CFD
%
% input:
% cfg_param: Cfg_filestr or class of 'SU2Config'.
% partitions: Processes number of parallel run SU2_CFD
% 
% output:
% SU2_data: Struct for ouput data of 'CONV_FILENAME'.csv file(if exist).
% SU2_history: Struct for history data of 'CONV_FILENAME'.csv file(if exist).
% SU2_surface: Struct for surface data of 'SURFACE_FILENAME'.csv file(if exist).
% SU2_CFD_info: Shell output of SU2_CFD.
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

% create SU2Config
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

if strcmp(config.getParameter('SOLVER'),'MULTIPHYSICS')
    error('Parallel computation script not compatible with MULTIPHYSICS solver.');
end

% get dir work
[str_filedir__,~]=fileparts(which('runSU2CFD.m'));
if isempty(dir_temp)
    dir_temp=fullfile(str_filedir__,'SU2_temp');
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
    error('runSU2CFD: config lack MESH_FILENAME define');
end
mesh_filestr=config.getParameter('MESH_FILENAME');
[mesh_filename,mesh_filedir,mesh_filestr]=safeSplitFileStr(mesh_filestr);
if ~exist(mesh_filestr,'file')
    error('runSU2CFD: mesh file do not exist')
end
config.setParameter('MESH_FILENAME',mesh_filename);
if contains(mesh_filename,'cgns','IgnoreCase',true)
    config.setParameter('MESH_FORMAT','CGNS');
elseif contains(mesh_filename,'su2','IgnoreCase',true)
    config.setParameter('MESH_FORMAT','SU2');
else
    error('runSU2CFD: unsupport mesh type');
end
safeCopyDirs(mesh_filename,mesh_filedir,dir_work,out_logger);

% process restart file
if config.isParameter('RESTART_FILENAME')
    restart_filestr=config.getParameter('RESTART_FILENAME');
    [restart_filename,restart_filedir,~]=safeSplitFileStr(restart_filestr);
    if exist(restart_filestr,'file')
        config.setParameter('RESTART_SOL','YES');
        config.setParameter('RESTART_FILENAME',restart_filename);
        safeCopyDirs(restart_filename,restart_filedir,dir_work,out_logger);
    end
end

% State
config.dump(fullfile(dir_work,'config_CFD.cfg'));
% state=SU2.io.State()
% out_logger.info(state)
if ~isempty(out_logger)
    out_logger.info(['mesh input file: ',mesh_filestr]);
    out_logger.info(['cfg file: ',cfg_filestr]);
    out_logger.info(['AOA: ',num2str(config.getParameter('AOA')),...
        ',SIDESLIP_ANGLE: ',num2str(config.getParameter('SIDESLIP_ANGLE')),...
        ',Ma: ',num2str(config.getParameter('MACH_NUMBER')),...
        ',T: ',num2str(config.getParameter('FREESTREAM_TEMPERATURE')),...
        ',P: ',num2str(config.getParameter('FREESTREAM_PRESSURE'))]);
end

% run SU2 CFD
% SU2.run.CFD(config)
if ~isempty(out_logger)
    out_logger.info('begin run SU2 CFD');
end

if ispc()
    run_command='SU2_CFD config_CFD.cfg';
else
    run_command=['mpirun -np ',num2str(partitions),' SU2_CFD config_CFD.cfg'];
end

% run command
cd(dir_work);
[status,SU2_CFD_info]=system(run_command);
cd(dir_cur);
if status ~= 0
    error_message=SU2_CFD_info;
else
    error_message=[];
    info_file=fopen(fullfile(dir_work,'SU2_CFD_info.log'),'w');
    fprintf(info_file,'%s',SU2_CFD_info);
    fclose(info_file);
    clear('info_file');
end

retry_time=0;
while (~isempty(error_message) && (retry_time < 2))
    if contains(error_message,'retrying')
        % retry again
        if ~isempty(out_logger)
            out_logger.warning('node busy,retry run SU2 CFD');
        end

        cd(dir_work);
        [status,SU2_CFD_info]=system(run_command);
        cd(dir_cur);
        if status ~= 0
            error_message=SU2_CFD_info;
        else
            error_message=[];
            info_file=fopen(fullfile(dir_work,'SU2_CFD_info.log'),'w');
            fprintf(info_file,'%s',SU2_CFD_info);
            fclose(info_file);
            clear('info_file');
        end

        retry_time=retry_time+1;
    else
        err_file=fopen(fullfile(dir_work,'SU2_CFD_err.log'),'w');
        fprintf(err_file,'%s',error_message);
        fclose(err_file);
        clear('err_file');
        if ~isempty(out_logger)
            out_logger.error([dir_work,error_message]);
        end
        error(['runSU2CFD: fatal error with SU2 CFD,proid=',num2str(procid)]);
    end
end

if ~isempty(out_logger)
    out_logger.info('end run SU2 CFD');
end

% obtain result
history_filestr=[config.getParameter('CONV_FILENAME'),'.csv'];
history_filestr=fullfile(dir_work,history_filestr);
if exist(history_filestr,"file")
    SU2_history=readSU2CSV(history_filestr);
    type_list=fieldnames(SU2_history);value_list=cell(length(type_list),1);
    for idx=1:length(value_list)
        value_list{idx}=SU2_history.(type_list{idx})(end);
    end
    SU2_data=cell2struct(value_list,type_list);
else
    SU2_history=[];
    SU2_data=[];
end

surface_filestr=[config.getParameter('SURFACE_FILENAME'),'.csv'];
surface_filestr=fullfile(dir_work,surface_filestr);
if exist(surface_filestr,"file")
    SU2_surface=readSU2CSV(surface_filestr);
else
    SU2_surface=[];
end

% delete temp file
if REMOVE_TEMP
    if ~isempty(out_logger)
        out_logger.info(['cleaning temp directory and files: ',dir_work]);
    end
    rmdir(dir_work,'s');
end

end
