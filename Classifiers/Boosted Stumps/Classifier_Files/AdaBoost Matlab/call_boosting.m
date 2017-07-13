function [trerr model] = callboosting(xtrain,ytrain,num_iter)
	
trerr=zeros(num_iter,1);
tsterr=zeros(num_iter,1);

model=boost_ada(xtrain,ytrain,num_iter);

for k=1:min(num_iter,size(model,1))
  y_trest=sign(eval_boost(model(1:k),xtrain));
  trerr(k) = sum(y_trest ~= ytrain);  
end

% plot(trerr./size(ytrain,1),'bo-')
% ylabel('Error');
% xlabel('Number of Boosting Iterations'), hold on,

% traccuracy = 1 - trerr(end)/size(ytrain,1)