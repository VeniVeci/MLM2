%% SCRIPT
% Manifold learning on the input data
% Export decomposition on filed input data.


clear

%% User Setup

nDimensionPreserve=1000;

manifoldLearning.name='PCA';

rawDataPath='./RawData/data_lambda_0_3_sigma_1_0/K_true.mat';
saveDataFolder='./PreprocessData/';

% f = fullfile(saveDataFolder,manifoldLearning.name,'kPCA.mat');

    
    %Addition setup for manifold learning
    switch manifoldLearning.name   
        case 'kPCA'
            manifoldLearning.ker='gaussian';
            manifoldLearning.new_dim = nDimensionPreserve;
    %         options.arg=10;   
            manifoldLearning.kAuto=1;

        case 'DiffusionMaps'
            manifoldLearning.metric ='euclidean';
            manifoldLearning.kernel ='gaussian'; 
            manifoldLearning.kAuto=1;
            manifoldLearning.dim_new = nDimensionPreserve-1;              
            manifoldLearning.t = 1;      
            manifoldLearning.Ztype = 1;
    %         options.FullRec = 0;      
        case 'Isomaps'
            manifoldLearning.dim_new=nDimensionPreserve;
        case 'PCA'
                       
        otherwise 
            error('No such DR method')
    end


%% Load data
load(rawDataPath);


%% manifold learning
switch manifoldLearning.name

    case 'kPCA'
        [Z,model] = Kpca(Ktrue',manifoldLearning);    
    case 'DiffusionMaps'
        [Z,model] = DiffusionMaps(Ktrue',manifoldLearning);   
    case 'Isomaps'
        [Z,model] = Isomaps(Ktrue',manifoldLearning);
    case 'PCA'
        [Z,model] = Lpca(Ktrue',nDimensionPreserve,manifoldLearning);
    otherwise 
        error('No such DR method')
end


%% 

SaveDataName=['K_' manifoldLearning.name '.mat'];
saveDataPath=fullfile(saveDataFolder,SaveDataName);
save(saveDataPath,'Z','model'); 


