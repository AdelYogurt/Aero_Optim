clc;
% clear;
close all hidden;

%% change all data

% for idx=1:1000
%     load([num2str(idx-1),'.mat']);
% 
%     str_type=[];
%     for str_idx=1:size(type,1)
%         str_type=[str_type,type(str_idx,:),','];
%     end
%     str_type=str_type(:)';
%     str_type=regexprep(str_type,'\s','');
%     str_type=regexprep(str_type,{'[',']','"'},'');
%     type_list=strsplit(str_type,',')';
% 
%     value_list=mat2cell(history,size(history,1),ones(1,size(history,2)))';
%     history=cell2struct(value_list,type_list(1:end-1));
% 
%     str_type=[];
%     for str_idx=1:size(surface.type,1)
%         str_type=[str_type,surface.type(str_idx,:),','];
%     end
%     str_type=str_type(:)';
%     str_type=regexprep(str_type,'\s','');
%     str_type=regexprep(str_type,{'[',']','"'},'');
%     type_list=strsplit(str_type,',')';
% 
%     value_list=mat2cell(surface.data,size(surface.data,1),ones(1,size(surface.data,2)))';
%     surface=cell2struct(value_list,type_list(1:end-1));
% 
%     save([num2str(idx-1),'.mat'],'AOA','Ma','history','surface','x')
% end


%% concentrate all data

% surface_list=[];
% value_list=[];
% x_list=[];
% for idx=1:1000
%     load([num2str(idx-1),'.mat']);
%     surface_list=[surface_list;surface];
%     value=history;
%     f=fieldnames(value);
%     for data_idx=1:length(f)
%         value.(f{data_idx})=value.(f{data_idx})(end);
%     end
%     value_list=[value_list;value];
%     x_list=[x_list;x];
% end
% 
% save('1000LHD.mat','AOA','Ma','value_list','surface_list','x_list')

%% generate database

load('1000LHD_light.mat','AOA','Ma','value_list','x_list')
CEff=zeros(1000,1);
Cl=zeros(1000,1);
X=x_list;
for idx=1:1000
    value=value_list(idx);
    Cl(idx)=value.CL;
    CEff(idx)=value.CEff;
end

up_bou=0.15*ones(1,14);
low_bou=[0.01,-0.05*ones(1,6),0.01,-0.05*ones(1,6)];

save('model_LDratio.mat','CEff','X','Cl','up_bou','low_bou');
