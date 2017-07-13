
function [h] = eval_stump(stump,X);

h = sign( stump.s*(X(:,stump.ind)-stump.x0) );
