function Coeff = pcasky(patterns_old)

%Reshape the data points using the principal component analysis
%
% Inputs:  patterns_old	- Input patterns
% Output:  Coeff: the eigenvector matix by PCA;
%
% Author: Mao Shasha (skymss0828@gmail.com),2008.6.15 %

if (nargin>3 || nargin<1) % check correct number of arguments
    help pcasky
else
    Coeff=[];
   if (nargin>=1 || nargin<=3)
       patterns=patterns_old';
       [r,c] = size(patterns);%%%注：r是特征数（维数）；c是样本个数；
       
       %Calculate cov matrix and the PCA matrixes
       m = mean(patterns')';
       S = ((patterns - m*ones(1,c)) * (patterns - m*ones(1,c))');%% S为协方差矩阵；     
       [Coeff latent]= eig(S);
   else
   end
end


