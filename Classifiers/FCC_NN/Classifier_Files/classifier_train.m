function [Decision_Surf]=classifier_train(data,group_no,varargin)
    [train_data, datamin, scale] = scaledata(data);
    targets = full(ind2vec(group_no'))';
    if nargin > 2
    hiddenLayerSize =  varargin{1};
    else
    hiddenLayerSize =3;
    end
    net =  cascadeforwardnet(hiddenLayerSize);
    net.trainFcn = 'trainbr';
    net.performFcn = 'mse';
    net.trainParam.epochs = 20;
    [net,tr] = train(net,train_data', targets');
    Decision_Surf = struct('net', net,'datamin', datamin,'scale_factor',scale);
