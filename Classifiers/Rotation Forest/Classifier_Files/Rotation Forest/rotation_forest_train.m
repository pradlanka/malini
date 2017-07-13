function [Decision_Surface]=rotation_forest_train(trainInp, trainY, K, ratio,L,knn)
numberfeature = size(trainInp,2);
classname= unique(trainY);
%%% training
%%% number of ensemble individuals;
for l=1:L
    %%% obtain the new samples by rotation forest %%%
    % K= floor(sqrt(numberfeature)); 
    % ratio=0.75;
    [R_new(:,:,l), R_coeff(:,:,l)]=RotationFal(trainInp, trainY, K, ratio);
    %%%% obtain new samples %%%%
    trainRFnew(:,:,l) = trainInp*R_new(:,:,l) ;  
end
Decision_Surface=struct('R_new',R_new, 'ratio', ratio, 'No_of_ensembles',L, 'Subsample_features',K, 'train_new',trainRFnew,'trainY',trainY,'knn',knn,'classname',classname);

