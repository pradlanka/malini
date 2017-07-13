function [testY]=rotation_forest_test(Decision_Surface, testInp)
L = Decision_Surface.No_of_ensembles;
trainRFnew = Decision_Surface.train_new;
trainY = Decision_Surface.trainY;
knn = Decision_Surface.knn;
R_new = Decision_Surface.R_new;
classname = Decision_Surface.classname;
%%% testing
for l=1:L
test_new = testInp*squeeze(R_new(:,:,l));
%%% learn for the new samples by classifiers(learners) %%%
%%% knn is the parameter of classifiers %%%
prelabeltest(:,l) = Nearest_Neighbor(squeeze(trainRFnew(:,:,l)), trainY, test_new,knn);
end
numbertest = size(testInp,1);
numberclass = size(classname,1);
%%% voting %%%
numberindex=[];
value=[];
testY=[];
for i=1:numbertest  
    prelabelES=[];
    prelabelES= prelabeltest(i,:); 
    for j=1:numberclass
        index=[];
        index=find(prelabelES==classname(j));
        numberindex(i,j)=length(index);
    end
    [value(i,1) indexmax(i,1)]=max(numberindex(i,:));
    testY(i,1)=classname(indexmax(i,1));
end