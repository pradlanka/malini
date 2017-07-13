function [performance]= runSVM_v3(data,groups,inner_kfold,itr)
%data: each arrow is one observation
%species : arrow(i) is the class member of observation data(i)
%clsname : The names of the g

% Select training and test sets randomly 
cp = classperf(groups);
for kk=1:itr
    indices = crossvalind('Kfold',length(groups),inner_kfold);
    for ij=1:inner_kfold
        test = (indices == ij); train =~ test;
        model = classifier_train(data(train,:),groups(train)');
        [classes] = classifier_predict(model,data(test,:));
        classperf(cp,classes,test);
    end
end
performance = cp.CorrectRate;
