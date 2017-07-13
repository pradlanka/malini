function [Decision_Surf]=classifier_train(data,group_no)
Decision_Surf = NaiveBayes.fit(data,group_no);