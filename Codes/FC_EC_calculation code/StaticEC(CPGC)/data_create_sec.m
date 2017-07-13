% If there are two classes Controls and PTSD they should be in the
% following format
addpath(genpath(pwd))
Filepath=pwd;
%X1=X{1,1};%staticFC_Controls;
%X2=X{1,2};%staticFC_PTSD;

coordinates=[]; count=0; 

cd(fullfile(folder_name,'output','subjects_together'));
Groupdir= dir;
y=cell(1,length(Groupdir)-2);
cd(Filepath);
for i=1:190
    for j= i+1:190
          count=count+1;
          for ipl =1:length(Groupdir)-2
          y{ipl} = cat(2,y{ipl},squeeze(X{1,ipl}(i,j,:)));
    end
          coordinates=[coordinates;[i,j]]; 
    end
end

for j=1:length(coordinates)
    connectivity_paths{j}= sprintf('''%d,%d''',coordinates(j,1),coordinates(j,2));
end
cd(fullfile(folder_name,'output','sec'))
Data =num2cell(vertcat(y{:,:}));
Table = horzcat(num2cell((1:length(coordinates))'), connectivity_paths',Data');

for psl=1:length(Groupdir)-2
    d=Groupdir(psl+2,1).name;xxx=size(X{1,psl},3);
    if (psl==1)
    labeletc=repmat({d},1,xxx);
    hore=labeletc;
    else
        labeletc=repmat({d},1,xxx);
        hore=horzcat(hore,labeletc);
    clear labeletc;
    end
    end
xlswrite('Table.xlsx', vertcat(num2cell((1:(size(Data,1)))),hore),'SFC','C1');
xlswrite('Table.xlsx',Table,'SFC','A3');

cd(Filepath);