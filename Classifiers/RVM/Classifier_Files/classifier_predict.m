function [classes_act] = classifier_predict(Decision_surf,data_test)
probs = Decision_surf.probs;
weights = Decision_surf.weights;
alphas = Decision_surf.alphas;
used = Decision_surf.used;
kernels_=Decision_surf.kernels;
P = Decision_surf.P;
kernel_width = Decision_surf.kernel_width;
kernel_data = Decision_surf.kernel_data;

PHI2=cell(P,1);
for p=1:P
    PHI2{p} = sbl_kernelFunction(data_test,kernel_data{p},kernels_{p},kernel_width);
end
predicted_vals = multinomial(PHI2,weights,size(data_test,1),P);
predicted = 1-sum(predicted_vals,2);
predicted=[predicted_vals, predicted];
[~, classes_act]= max(predicted,[],2);
