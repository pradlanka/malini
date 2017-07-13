function variable_importance = find_predictorImportance(Decision_surface,spot_id_tmp,spotid)
variable_importance =  zeros(length(spotid),1);
selected_paths = ismember(spotid,spot_id_tmp);  %Point to the paths been selected to the analysis
importance_scores = predictorImportance(Decision_surface);
variable_importance(selected_paths)= importance_scores;
