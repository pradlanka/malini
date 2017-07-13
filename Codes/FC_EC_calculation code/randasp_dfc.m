Filepath=pwd;
cd(fullfile(folder_name,'output','dfc'));
[~,~,ltinput_data]=xlsread('Table.xlsx','SFC');
%YourFile=handles.YourFile;
%ltinput_data=YourFile;

%classname{2}='EMCI';
%classname{3}='LMCI';
%classname{4}='AD';
%grp=classname;
nocl=handles.nocl;
no_classes=nocl;
[row, col] = size(ltinput_data);         
data = cell2mat(ltinput_data(3:row,3:col));
samples_name  = cell2mat(ltinput_data(1,3:col));                % Extract subjects 
[~, n_samples]= size(samples_name);
groups_name   = ltinput_data(2,3:col);% Extract the group names
classname=unique(groups_name);
grp=classname;
paths = ltinput_data(3:row,2);
spot_id = cell2mat(ltinput_data(3:row,1));
classes = 1:length(unique(groups_name));
% Convert Group names into Group Indices
groups_ID=zeros(n_samples,1);
for ns=1:n_samples
    for ng=1:no_classes
        if(strcmp(groups_name(ns),classname{ng})==1)
        groups_ID(ns)=ng;
        end
    end
end
[grpnam,loc]=unique(groups_ID);
edc=ltinput_data(:,1);
cdex=ltinput_data(:,2);
onedata=ltinput_data(:,3:end);
for k=1:no_classes
    if(k==no_classes)
        grpdat{k}=onedata(:,loc(k):end);
    else
grpdat{k}=onedata(:,loc(k):loc(k+1)-1);
    end
end
splirat=handles.splirat;num_points = size(onedata,2);
split_point = round(num_points*splirat);
for l=1:no_classes
    num_points = size(grpdat{l},2);
split_point = round(num_points*splirat);
    grpdats=grpdat(l);
    grpdatind=grpdats{1,1};
X_train{l} = grpdatind(:,1:split_point);
X_test{l} =grpdatind(:,split_point+1:end);
end
hope=cat(2,edc,cdex,X_train{1,1});
hopet=cat(2,edc,cdex,X_test{1,1});
%inputhope_data=cat(2,hope,X_train{1,s});
%testhope_data=cat(2,hopet,X_test{1,s});
for s=2:no_classes
hope=cat(2,hope,X_train{1,s});
hopet=cat(2,hopet,X_test{1,s});
end
input_data=hope;
test_data=hopet;
handles.input_data=input_data;
assignin('base','input_data',handles.input_data);
guidata(hObject,handles);
handles.test_data=test_data;
assignin('base','test_data',handles.test_data);
guidata(hObject,handles);
cd(Filepath);