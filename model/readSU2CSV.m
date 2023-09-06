function data=readSU2CSV(filestr)
% read analysis data of SU2
%
if ~exist(filestr,'file')
    data=[];
    return;
end

% read out type
file_data=fopen(filestr,"r");
str_type=regexprep(fgetl(file_data),'\s','');
str_type=regexprep(str_type,{'[',']','"','-'},'');
type_list=strsplit(str_type,',')';
fclose(file_data);
clear('file_data');

% read value
value_list=importdata(filestr,',',1).data;
value_list=mat2cell(value_list,size(value_list,1),ones(1,size(value_list,2)))';

data=cell2struct(value_list,type_list);
end
