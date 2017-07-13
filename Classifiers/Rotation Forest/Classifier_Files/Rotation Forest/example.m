%%%%  an example for rotation forest %%%%
clc;
clear all;

load blfdna.mat
trainXOrg=dnasamples;
trainYOrg=dnalabels;
testXOrg=dnatestsamples;
testYOrg=dnatestlabels;
numberfeature=size(trainXOrg,2);
numbertrain=length(trainYOrg);
numbertest=length(testYOrg);

% load blfwdbc.mat
% dataX=wdbcsamples;
% dataY=wdbclabels;
% [number numberfeature]=size(dataX);
% numbertrain=100;
% numbertest=number-numbertrain;
% trainXOrg=dataX(1:numbertrain,:);
% trainYOrg=dataY(1:numbertrain,:);
% testXOrg=dataX(1+numbertrain:end,:);
% testYOrg=dataY(1+numbertrain:end,:);
%%% number of classes of samples %%
class=unique(testYOrg);
numberclass=length(class);
%%% training
[Decision_Surface]=rotation_forest_train(trainXOrg, trainYOrg, 13, 0.75, 6,13);
[preENCRF]=rotation_forest_test(testXOrg,Decision_Surface,class);
%%% compute the accuracy rate of ensemble %%%
accuracyrate=sum(preENCRF==testYOrg)/numbertest;
        