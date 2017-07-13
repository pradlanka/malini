function [classes_act] = classifier_predict(Decision_surf,data_test)
[classes_act] = elm_kernel_predict(Decision_surf,data_test);