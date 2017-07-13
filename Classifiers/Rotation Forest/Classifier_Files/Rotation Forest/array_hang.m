
function [data_new_hang index_hang]=array_hang(X)
%%% Author: Mao Shasha,skymss0828@gmail.com,2008.6.2 %%%%

data_ini=[];
data_ini=X;
index_hang=[];
data_new_hang=[];

%%% for row vector of X %%%
rand('state',sum(100*clock)); 
index_hang=randperm(size(data_ini,1));    
data_new_hang=data_ini(index_hang,:);
