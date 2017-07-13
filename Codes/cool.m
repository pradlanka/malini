addpath(genpath(pwd));
clas2 = get(handles.elmcheck,'Value');
clas3 = get(handles.rceknncheck,'Value');
clas4 = get(handles.ldacheck,'Value');
clas5 = get(handles.svmcheck,'Value');
clas6 = get(handles.naivecheck,'Value');
clas7 = get(handles.qdacheck,'Value');
clas8 = get(handles.rbfcheck,'Value');
clas9 = get(handles.baggedcheck,'Value');
clas10 = get(handles.stumpscheck,'Value');
clas11 = get(handles.boostedcheck,'Value');
clas12 = get(handles.fccnncheck,'Value');
%clas13 = get(handles.oknncheck,'Value');
clas14 = get(handles.lvqnetcheck,'Value');
clas15 = get(handles.mlpcheck,'Value');
clas16 = get(handles.randomcheck,'Value');
clas17 = get(handles.rlrcheck,'Value');
clas18 = get(handles.rotationcheck,'Value');
clas19 = get(handles.rvmcheck,'Value');
clas20 = get(handles.slrcheck,'Value');
clas21 = get(handles.consencheck,'Value');
rceelm = get(handles.checkbox55,'Value');
rceknn = get(handles.checkbox57,'Value');
rcelda = get(handles.checkbox58,'Value');
rcelinearsvm = get(handles.checkbox59,'Value');
rcenaive = get(handles.checkbox60,'Value');
rceqda = get(handles.checkbox61,'Value');
rcebagged = get(handles.checkbox68,'Value');
rceboosted = get(handles.checkbox67,'Value');
rcetrees = get(handles.checkbox66,'Value');
rcefccnn = get(handles.checkbox65,'Value');
%rcewknn = get(handles.checkbox64,'Value');
rcemlp = get(handles.checkbox63,'Value');
rceslr = get(handles.checkbox62,'Value');
rceratfo = get(handles.checkbox69,'Value');
rcerbf = get(handles.checkbox70,'Value');
rcerlr = get(handles.checkbox71,'Value');
rcelvq = get(handles.checkbox72,'Value');
rcerandfo = get(handles.checkbox73,'Value');
rcervm = get(handles.checkbox74,'Value');
%rceelm = get(handles.checkbox55,'Value');

loaddata = get(handles.checkbox34,'Value');
if(loaddata==1)
random_split_proper;
end

if(clas2==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
    discos=strcat(idod,'ELM');
    arerce=rceelm;
    main_rce_elm;
    Test_surface_elm;
    clearvars arerce discos;
    end
    if(clas3==1)
         global discos;
         if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'KNN');
     arerce=rceknn;
    main_rce_knn;
    Test_surface_rce_knn;
    clearvars discos arerce;
    end
    if(clas4==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'LDA');
     arerce=rcelda;
    main_lda;
    Test_surface_lda;
    clearvars discos arerce;
    end
    if(clas5==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
        discos=strcat(idod,'Linear_SVM');
        arerce=rcelinearsvm;
 main_linear_svm;
    Test_surface_linear_svm;
    clearvars discos arerce;
    end
    if(clas6==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'NAIVE_BAYES');
     arerce=rcenaive;
    main_naivebayes;
      Test_surface_naive;
    clearvars arerce discos;
    end
    if(clas7==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'QDA');
     arerce=rceqda;
        main_qda;
    Test_surface_qda;
    clearvars discos arerce;
    end
    if(clas8==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
        discos=strcat(idod,'RBF_SVM');
    arerce=rcerbf;
        main_rbf_svm;
    Test_surface_svm;
    clearvars arerce discos;
    end
    if(clas9==1)
         global discos;
         if(loaddata==1)
            idod='independent';
        end
        discos=strcat(idod,'Bagged_Trees');
      arerce=rcebagged;
       main_bagged;
    Test_surface_bagged
    clearvars discos arerce;
       end
    if(clas10==1)
         global discos;
     if(loaddata==1)
            idod='independent';
        end
        discos=strcat(idod,'Boosted_stumps');
    arerce=rceboosted;
       main_boosted_stumps;
    Test_surface_boosted_stumps;
    clearvars discos arerce;
    end
    if(clas11==1)
         global discos;
         if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'BOOSTED_TREES');
    arerce=rcetrees;
     main_boosted;
    Test_surface_boosted_trees;
    clearvars discos arerce;
    end
    if(clas12==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'FCC_NN');
    arerce=rcefccnn;
     main_fcc_nn;
    Test_surface_fcc_nn;
    clearvars discos arerce;
    end
   % if(clas13==1)
    %    global discos;
     %    if(loaddata==1)
      %      idod='independent';
       % end
        %discos=strcat(idod,'KNN');
    %arerce=rcewknn;
     %   main_knn;
    %Test_surface_knn;
    %clearvars discos arerce;
    %end
    if(clas14==1)
         global discos;
         if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'LVQNET');
    arerce=rcelvq;
     main_lvqnet;
    Test_surface_lvqnet;
    clearvars discos arerce;
    end
    if(clas15==1)
       global discos;
       if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'MLP');
     arerce=rcemlp;
    main_mlp;
    Test_surface_mlp;
    clearvars discos arerce;
    end
    if(clas16==1)
         global discos;
         if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'RANDOM_FOREST');
    arerce=rcerandfo;
     main_random_forest;
    Test_surface_rando;
    clearvars discos arerce;
    end
    if(clas17==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'RLR');
    arerce=rcerlr;
     main_rlr;
    Test_surface_rlr;
    clearvars discos arerce;
    end
    if(clas18==1)
        global discos;
        if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'ROTATION_FOREST');
    arerce=rceratfo;
     main_rotfo;
    Test_surface_rotfo;
    clearvars discos arerce;
    end
    if(clas19==1)
         global discos;
         if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'RVM');
    arerce=rcervm;
     main_rvm;
    Test_surface_rvm;
    clearvars discos arerce;
    end
    if(clas20==1)
          global discos;
          if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'SLR');
     arerce=rceslr;   
     main_slr;
        Test_surface_slr;
        clearvars arerce discos;
    end
    if(clas21==1)
          global discos;
          if(loaddata==1)
            idod='independent';
        end
     discos=strcat(idod,'Combined_Results');
       CombinedResults_new;
        clearvars discos;
    end
     