%load test_data.mat
currentfolder=pwd;
%[~,~,test_data]=xlsread('train_data.xlsx');
test_data=handles.test_data;
test_data_targets = (test_data(2,3:end)); 
coordinates = test_data(3:end,2);
test_data=(test_data(3:end,3:end));
test_data_inp =(cell2mat(test_data))';
%classname{1}='Controls';
%classname{2}='EMCI';
% Convert group names to group nos for the testing data 
no_classes=size(classname,2);  
group_no_targets=zeros(size(test_data_targets)) ;
for ns=1:1:length(test_data_targets);
    for ng=1:no_classes
        if(strcmp(test_data_targets(ns),classname{ng})==1)
        group_no_targets(ns)=ng;
        end
    end
end
 cp2 = classperf(group_no_targets);
classes  = 1:no_classes;
cd(fullfile(folder_name,'Outside_Rce_Results'));  % CLassifiers implented outside the RCE framework
Algs_group1 = dir; 
Test_Results_table1 = []; bal_acc_classifier1 = [];
for  classifiers1 = 3:length(Algs_group1)
   cd(fullfile(folder_name,'Outside_Rce_Results'));
if(Algs_group1(classifiers1).isdir)
cd(Algs_group1(classifiers1).name)

% Load a few variables from  Results
%cd(fullfile(currentfolder,'results_bagged'));
load Test_results.mat
%classname{1}='Controls';
%classname{2}='EMCI';
classprobs = zeros(no_classes,length(test_data_targets));
classaccuracy = zeros(no_classes,1);
for j=1:length(test_data_targets)
    targetvals= squeeze(classes_predicted(:,j));
    for n=1:no_classes
        classprobs(n,j) = sum(targetvals==n);
    end
     classprobs(:,j) =  classprobs(:,j)/sum(classprobs(:,j));
end
for j=1:no_classes
    classaccuracy(j) = Conf_mat_train(j,j)/sum(squeeze(Conf_mat_train(j,:)));
end
bal_accuracy =mean(classaccuracy);
Test_Results_table1 = cat(3,Test_Results_table1,classprobs);
 bal_acc_classifier1 = [bal_acc_classifier1 bal_accuracy];
cd ..
end
end
cd ..
cd(fullfile(folder_name,'Rce_Results'));  % For classifiers implented outside the RCE framework
Algs_group2 = dir; 
Test_Results_table2 = []; bal_acc_classifier2=[];
for  classifiers2 = 3:length(Algs_group2)
    cd(fullfile(folder_name,'Rce_Results'));
if(Algs_group2(classifiers2).isdir)
cd(Algs_group2(classifiers2).name)

% Load a few variables from  Results
%cd(fullfile(currentfolder,'rce_results'));
load Test_results.mat
classprobs = zeros(no_classes,length(test_data_targets));
classaccuracy = zeros(no_classes,1);
for j=1:length(test_data_targets)
    targetvals= squeeze(classes_predicted(:,j));
    for n=1:no_classes
        classprobs(n,j) = sum(targetvals==n);
    end
end
for j=1:no_classes
    classaccuracy(j) = Conf_mat_train(j,j)/sum(squeeze(Conf_mat_train(j,:)));
end
 bal_accuracy =mean(classaccuracy);
 Test_Results_table2 = cat(3,Test_Results_table2,classprobs);
 bal_acc_classifier2 = [bal_acc_classifier2 bal_accuracy];
%cd ..
%cd ..
end
end
cd ..
Test_Results_AllClassifiers  =cat(3,Test_Results_table1,Test_Results_table2);
 bal_acc_classifier_all = [bal_acc_classifier1 bal_acc_classifier2];
clearvars -except Test_Results_AllClassifiers group_no_targets classes classname cp2 no_classes bal_acc_classifier_all
Best_Results =zeros(length(group_no_targets),1);
for nl= 1:length(group_no_targets)
    [~,Best_Results(nl)] = max(squeeze(Test_Results_AllClassifiers(:,nl,:))*bal_acc_classifier_all');
end
classperf(cp2,Best_Results);
Confmat_test=confusionmat(group_no_targets,Best_Results,'order',classes);
for j=1:no_classes
    assignin('base',strcat('acc_test_',classname{j}), Confmat_test(j,j)/sum(squeeze(Confmat_test(j,:))));
end
global folder_name;
cd(folder_name);
mkdir('Consesus_results');
cd('Consesus_results');
Combined_Best_accuracy = cp2.CorrectRate;
save('Consensus_Results.mat');
cd(currentfolder);
