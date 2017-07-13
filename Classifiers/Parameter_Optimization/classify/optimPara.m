function [performance, specificity, sensitivity, Conf_mat, Decision_Surf_best, bestparam]= optimPara(data,group_no,data_test, group_no_test,kfold,resample,classes,Hyperparameter)
% data: each arrow is one observation

cp2 = classperf(group_no_test);

% Search for the best set of parameters
Optimizable_parameters =length(Hyperparameter);
Accuracy_paramset = zeros(prod(sizeout(Hyperparameter)),1);
Paramlist = paramgrid(Hyperparameter);
for i=1:prod(sizeout(Hyperparameter))
    Paramset =Paramlist(i,:);
    cp = classperf(group_no);
    for j=1:resample
    indices = crossvalind('Kfold',group_no,kfold);
    for ij=1:kfold
        test = (indices == ij); train =~ test;
        Decision_Surf = classifier_train(data(train,:),group_no(train),Paramset);
        [classes_act] = classifier_predict(Decision_Surf, data(test,:));
        classperf(cp,classes_act,test);
    end
    end
    Accuracy_paramset(i)=cp.CorrectRate;
end
[max_Accuracy, Index] =  max(Accuracy_paramset);
bestparam= Paramlist(Index,:);

% Train using the best value of parameters
Decision_Surf_best = classifier_train(data, group_no, bestparam);
[classes_act_test] = classifier_predict(Decision_Surf_best,data_test);

% Updates the CP object with the current classification results
classperf(cp2,classes_act_test);
performance = cp2.CorrectRate;
specificity = cp2.specificity;
sensitivity = cp2.sensitivity;

% Calculate performance Metrics
Conf_mat=confusionmat(group_no_test,classes_act_test,'order',classes);


function [no_para] = sizeout(x)
nout = length(x);no_para =zeros(1,nout);
for k = 1:nout
   no_para(k) = length(x{k});
end