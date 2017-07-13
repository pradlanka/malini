function X=CreateData

X1=[0.6*rand(50,1) 0.7*rand(50,1) 1*ones(50,1)];
X2=[1.0-0.6*rand(50,1) 0.6*rand(50,1) 2*ones(50,1)];
X3=[1.0-0.7*rand(50,1) 1.0-0.6*rand(50,1) 3*ones(50,1)];
X=[X1;X2;X3];


