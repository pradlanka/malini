
%%% skymss,2008.6.15 (email: skymss@126.com) %%%%

function [R_new,R_coeff]=RotationFal(Xini, Yini, K, Ratio)

% RotationFal: obtaining the new samples by Rotation forest algorithm;
%
% function [R_new,R_coeff]=RotationFal(DataX_old,trainY,K,Ratio)
%
% Input:  Xini: the original samples;
%         Yini: the labels of original samples;
%            K: the number of feature subset of Xini;
%         Ratio: the proportion of resampling new samples from original
%         samples; (usally Ratio=75%)
% Ourput: R_new: the arranged coffecients of features by Rotation forest;
%         R_coeff: the unarranged coffecients of features by Rotation forest;
%
% Author: Mao Shasha (skymss0828@gmail.com),2008.6.15 %


if (nargin<4 || nargin>5)
    help RotationFal
else
    Xinirfnew=[];

    DataX=[];
    DataX=Xini;
    [numtrn num_feature]=size(DataX);
    M=floor(num_feature/K);
    class=[];
    class=unique(Yini);
   
    [data_new_lie index_lie_new]=array_lie(DataX);
    trainX_lie=[];
    trainX_lie=data_new_lie;
    R_coeff=zeros(num_feature,num_feature);
    trainX_subset=[];
    trainY=Yini;
    index=[];
    for m=1:length(class)
        index{m,1}=find(trainY==class(m));
    end
    %%%%%%%%%%%%%%%%%%%%%%% 划分特征样本进行PCA转换 %%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:K
        number1=(i-1)*M;
        number2=i*M;
        if (number2<=num_feature)
            trainX_subset=trainX_lie(:,number1+1:number2);
        else
            trainX_subset=trainX_lie(:,number1+1:end);
        end
        %%%% 抽取非空的类别子集 %%%%
        %%%% eliminate from dataset a random subset of classes %%%
        rateeliminate=0.08;
        [AsubX,AsubY]=randomsub([trainX_subset,trainY],rateeliminate);
        %%%% using bootstrap algorithm to obtain subset of samples (进行bootstrap抽样) %%%%
        [trainX_subset_new,trainY_new,indexselect]=bootstrapal(AsubX,AsubY,Ratio);
        %%%% using PCA to transform samples (进行PCA处理，产生新样本(注：这里的PCA处理样本是特征*样本数)) %%%%
        Coeff=pcasky(trainX_subset_new); 
        A_coeff{i,1}=Coeff; 
        number3=number1+size(Coeff,2);
        R_coeff(number1+1:number2,number1+1:number3)=Coeff;
    end
    %%%% arrage R_coeff based on original feature (对R_coeff重新排列获得R_new) %%%%
    [index_lie index_A]=sort(index_lie_new);
    R_new1=R_coeff(:,index_A);
    R_new=[];
    for i=1:num_feature
        number_zero=length(find(R_new1(:,i)==0));
        if (number_zero>=size(R_new1,1))
            R_new=R_new;
        else
            R_new=[R_new R_new1(:,i)];    
        end
    end

    %%%% obtain new samples %%%%
    Xinirfnew=Xini*R_new;       
end

