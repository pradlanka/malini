function [classes_act] = classifier_predict(Decision_surf,data_test)
net = Decision_surf.net;
datamin = Decision_surf.datamin;
scale = Decision_surf.scale_factor;
test_data = scaledata(data_test',datamin,scale);
tstOutputs = net(test_data);
[~,Index]=max(tstOutputs);
classes_act= Index';