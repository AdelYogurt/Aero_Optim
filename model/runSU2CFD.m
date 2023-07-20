function [SU2_history,SU2_surface,SU2_CFD_info]=runSU2CFD...
    (mesh_filestr,cfg_filestr,partitions,CFD_config,restart_filestr,...
    dir_temp,run_description,remove_dump,out_logger)
% interface of SU2_CFD
% base on input filestr and CFD_config to run SU2_CFD
%
if nargin < 9
    out_logger=[];
    if nargin < 8
        remove_dump=[];
        if nargin < 7
            run_description=[];
            if nargin < 6
                dir_temp=[];
                if nargin < 5
                    restart_filestr=[];
                    if nargin < 4
                        CFD_config=[];
                    end
                end
            end
        end
    end
end

[str_filedir__,~]=fileparts(which('runSU2CFD.m'));

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
if ~isempty(restart_filestr)
    [restart_filename,restart_filedir,restart_filestr]=safeSplitFileStr(restart_filestr);
end

% Config
config=Config(cfg_filestr);
config.setParameter('NUMBER_PART',partitions);

if strcmp(config.getParameter('SOLVER'),'MULTIPHYSICS')
    error('Parallel computation script not compatible with MULTIPHYSICS solver.');
end

% process input
if ~isempty(restart_filestr)
    config.setParameter('RESTART_SOL','YES');
    config.setParameter('RESTART_FILENAME',restart_filename);
end
config.setParameter('MESH_FILENAME',mesh_filename);
if contains(mesh_filename,'cgns')
    config.setParameter('MESH_FORMAT','CGNS');
else
    config.setParameter('MESH_FORMAT','SU2');
end

% add input config
if ~isempty(CFD_config)
    field_name_list=fieldnames(CFD_config);
    for idx=1:length(field_name_list)
        field=field_name_list{idx};
        config.setParameter(field,CFD_config.(field));
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
if ~isempty(restart_filestr)
    % copy restart data to dir work
    safeCopyDirs({restart_filename},restart_filedir,dir_work,out_logger);
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
    run_command=['mpirun -np ',num2str(config.getParameter('NUMBER_PART')),' SU2_CFD config_CFD.cfg'];
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
            out_logger.error([dir_work,error_message,'\n']);
        end
        error(['runSU2CFD: fatal error with SU2 CFD,proid=',num2str(procid)]);
    end
end

if ~isempty(out_logger)
    out_logger.info('end run SU2 CFD');
end

% obtain result
history_filestr=[config.getParameter('CONV_FILENAME'),'.csv'];
SU2_history=readSU2CSV(fullfile(dir_work,history_filestr));
surface_filestr=[config.getParameter('SURFACE_FILENAME'),'.csv'];
SU2_surface=readSU2CSV(fullfile(dir_work,surface_filestr));

% delete temp file
if remove_dump
    if ~isempty(out_logger)
        out_logger.info('cleaning temp files...');
    end
    rmdir(dir_work,'s');
end

end
