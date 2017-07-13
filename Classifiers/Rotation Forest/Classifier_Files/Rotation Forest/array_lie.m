

function [data_new_lie index_lie]=array_lie(X)
%%% Author: Mao Shasha,skymss0828@gmail.com,2008.6.2 %%%%

data_ini=[];
data_ini=X;
data_new_lie=[];
index_lie=[];

%%% for column vector of X %%%%
rand('state',sum(100*clock)); 
index_lie=randperm(size(data_ini,2));    
data_new_lie=data_ini(:,index_lie);
