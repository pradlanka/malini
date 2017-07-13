function [performance, specificity, sensitivity, Conf_mat, Decision_Surf]= runclassifier_Default(data,group_no,data_test, group_no_test,classes,Decision_Surface_Flag)

 cp = classperf(group_no_test);

% Train using the best value of parameters
Decision_Surf = classifier_train(data, group_no);
[classes_act_test] = classifier_predict(Decision_Surf,data_test);

% Updates the CP object with the current classification results
classperf(cp, classes_act_test);
performance = cp.CorrectRate;
specificity = cp.specificity;
sensitivity = cp.sensitivity;

% Calculate performance Metrics
Conf_mat=confusionmat(group_no_test,classes_act_test,'order',classes);
