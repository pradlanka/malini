function [Decision_Surf]=classifier_train(data,group_no)
Decision_Surf = fitcdiscr(data, group_no,'DiscrimType','quadratic','Gamma',1);