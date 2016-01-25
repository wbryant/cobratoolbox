function [AddedRxns] = submitFastGapFill(inputsFile,modelFileIn,dbFileIn,dictionaryFileIn,workspaceFileIn,forceRerun)
%% function [AddedRxns] = submitFastGapFill(modelFile,dbFile,dictionaryFile,workspaceFile,paramsFile)
% William Bryant, Jan 2016
%


%% Preparation - load data

% Suppress load warnings and deal with in the script
warning('off','MATLAB:load:variableNotFound')

% Default files
modelFile = '/Users/wbryant/work/BTH/analysis/fastGapFill/defaults/input/BTH_iAH991_w_gprs.xml';
dbFile = '/Users/wbryant/work/BTH/analysis/fastGapFill/defaults/input/BTH_iAH991_w_gprs.xml';
dictionaryFile = '/Users/wbryant/work/BTH/analysis/fastGapFill/defaults/input/BTH_iAH991_w_gprs.xml';
workspaceFile = '/Users/wbryant/work/BTH/analysis/fastGapFill/defaults/input/BTH_iAH991_w_gprs.xml';

if exist('inputsFile','var') && ~isempty('inputsFile')
    %Get input filenames from inputsFile
    try
        fprintf('Obtaining parameters from inputsFile ...\n')
        fileHandle = fopen(inputsFile);
        u = textscan(fileHandle,'%s\t%s');
        for i = 1:length(u{1})
            paramName = matlab.lang.makeValidName(u{1}{i});
            paramValue = u{2}{i};
            eval([paramName ' = ''' paramValue ''';'])
        end
        fclose(fileHandle);
    catch
        fprintf('Parameters could not be read from file, using defaults ...\n')
    end
else
    % Get input files where specified
    if exist('modelFileIn','var') && ~isempty(modelFileIn)
        modelFile = modelFileIn;
    end
    if exist('dbFileIn','var') && ~isempty(dbFileIn)
        dbFile = dbFileIn;
    end
    if exist('dictionaryFileIn','var') && ~isempty(dictionaryFileIn)
        dictionaryFile = dictionaryFileIn;
    end
    if exist('workspaceFileIn','var') && ~isempty(workspaceFileIn)
        workspaceFile = workspaceFileIn;
    end
end


% baseDirectory = '/Users/wbryant/work/BTH/analysis/fastGapFill/defaults/';
% inputDir = 'input/';
% modelFile = strcat(baseDirectory,inputDir,'BTH_iAH991_w_gprs.xml');
% dbFile = strcat(baseDirectory,inputDir,'reaction.lst');
% dictionaryFile = strcat(baseDirectory,inputDir,'KEGG_dictionary.xls');
% workspaceFile = strcat(baseDirectory,'output/','default_out.mat');




forceRerun=true;
   
oldDirectory = cd(baseDirectory);
epsilon = [];

% If workspaceFile is present, check for all variables; if consistModel, 
% consistMatricesSUX and BlockedRxns are present, do not rerun 
% prepareFastGapFill, otherwise run and save workspace; check for but do
% not require prepStats
try
    a = load(workspaceFile,'consistModel','consistMatricesSUX','BlockedRxns','model','baseDirectory');
catch
    fprintf('Workspace file not found, proceeding with prepareGapFill ...\n')
    a = struct();
end
length_a = length(fieldnames(a));
clear a;
if (length_a == 5) && ~forceRerun
    fprintf('prepareGapFill already run, proceed with run.\n');
else
% load model using relevant load function
    if regexp(modelFile,'.mat$')
        model = readMlModel(modelFile);
    elseif regexp(modelFile,'.xml$')
        model = readCbModel(modelFile);
    end
    
    % Switch all IDs to lower-case
    model.mets = cellfun(@lower,model.mets,'UniformOutput',false);

    % remove constraints from exchange reactions
    EX = strncmp('EX_',model.rxns,3);
    model.lb(EX)=-100;
    model.ub(EX)=100;
    clear EX

    tic;
    [consistModel,consistMatricesSUX,BlockedRxns] = prepareFastGapFill(model, {}, epsilon, dbFile, dictionaryFile);
    tpre=toc;

    % Prepare the output table with statistics
    cnt = 1;
    prepStats{cnt,1} = 'Model name';cnt = cnt+1;
    prepStats{cnt,1} = 'Size S (original model)';cnt = cnt+1;
    prepStats{cnt,1} = 'Number of compartments';cnt = cnt+1;
    prepStats{cnt,1} = 'List of compartments';cnt = cnt+1;
    prepStats{cnt,1} = 'Number of blocked reactions';cnt = cnt+1;
    prepStats{cnt,1} = 'Number of solvable blocked reactions';cnt = cnt+1;
    prepStats{cnt,1} = 'Size S (flux consistent)';cnt = cnt+1;
    prepStats{cnt,1} = 'Size SUX (including solvable blocked reactions)';cnt = cnt+1;
    prepStats{cnt,1} = 'Time preprocessing';

    % get stats
    cnt = 1;
    prepStats{cnt,2} = modelFile;cnt = cnt+1;
    [a,b] = size(model.S);
    prepStats{cnt,2} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
    [~,rem] = strtok(model.mets,'\[');
    rem = unique(rem);
    prepStats{cnt,2} = num2str(length(rem));cnt = cnt+1;
    Rem = rem{1};
    for j = 2:length(rem)
        Rem = strcat(Rem,',',rem{j});
    end
    prepStats{cnt,2} = Rem;cnt = cnt+1;
    clear Rem rem;
    prepStats{cnt,2} = num2str(length(BlockedRxns.allRxns));cnt = cnt+1;
    prepStats{cnt,2} = num2str(length(BlockedRxns.solvableRxns));cnt = cnt+1;
    [a,b] = size(consistModel.S);
    prepStats{cnt,2} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
    [a,b] = size(consistMatricesSUX.S);
    prepStats{cnt,2} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
    prepStats{cnt,2} = tpre;cnt = cnt+1;
    
    % Save data then clear all except required variables
    save(workspaceFile,'consistModel','consistMatricesSUX','BlockedRxns','prepStats','model','baseDirectory');
    cd(oldDirectory);
    clear
end

fprintf('Running fastGapFill ...\n')
% define weights for reactions to be added - the lower the weight the
% higher the priority


%% Define runs
% WARNING: weights should not be 1, as this impacts on the algorithm
% performance

runs = {};

%% RUN 1 
run.weightsPerReactionFile = 'metabolonWeightedRxns.tsv';
run.weights.MetabolicRxns = 0.1; % Kegg metabolic reactions
run.weights.ExchangeRxns = 0.5; % Exchange reactions
run.weights.TransportRxns = 10; % Transport reactions
run.name = 'initial';
runs = [runs; run];

%% RUN 2
run.weightsPerReactionFile = 'metabolonWeightedRxns.tsv';
run.weights.MetabolicRxns = 0.1; % Kegg metabolic reactions
run.weights.ExchangeRxns = 0.9; % Exchange reactions
run.weights.TransportRxns = 40; % Transport reactions
run.name = 'high_weight';
runs = [runs; run];

clear run

%% Do not change below here
% Prepare the output table with statistics
% Stats{cnt,1} = 'Model name';cnt = cnt+1;
% Stats{cnt,1} = 'Size S (original model)';cnt = cnt+1;
% Stats{cnt,1} = 'Number of compartments';cnt = cnt+1;
% Stats{cnt,1} = 'List of compartments';cnt = cnt+1;
% Stats{cnt,1} = 'Number of blocked reactions';cnt = cnt+1;
% Stats{cnt,1} = 'Number of solvable blocked reactions';cnt = cnt+1;
% Stats{cnt,1} = 'Size S (flux consistent)';cnt = cnt+1;
% Stats{cnt,1} = 'Size SUX (including solvable blocked reactions)';cnt = cnt+1;

cnt = 1;
clear Stats;
Stats{cnt,1} = 'Run name';cnt = cnt+1;
Stats{cnt,1} = 'Number of added reactions (all)';cnt = cnt+1;
Stats{cnt,1} = 'Number of added metabolic reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added transport reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added exchange reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Time fastGapFill';cnt = cnt+1;
Stats{cnt,1} = 'SFilename';cnt = cnt+1;

workspaceFile = '/Users/wbryant/work/BTH/analysis/fastGapFill/iAH991/output/iAHPrepared.mat';
a = load(workspaceFile,'consistModel','consistMatricesSUX','BlockedRxns','model','baseDirectory');
length_a = length(fieldnames(a));
clear a;
if length_a == 5
    load(workspaceFile,'consistModel','consistMatricesSUX','BlockedRxns','model','baseDirectory');
    fprintf('Variables loaded from .mat file ...\n');
    oldDirectory = cd(baseDirectory);
    b = load(workspaceFile,'prepStats');
    if isempty(b)
        load(workspaceFile,'prepStats');
        fprintf('Warning: no stats for prepareFastGapfill ...');
    end
    clear b
else
    fprintf('prepareFGF must be run first ...')
    return
end


baseDirectory=strcat(baseDirectory,'/');

for i = 1:length(runs)
    
     
    run_idx = i+1;
    fprintf('\nRun %i (%s)\n',i,runs{i}.name);
    weights = runs{i}.weights;
    weightsPerReactionFile = runs{i}.weightsPerReactionFile;
    SFilename = strcat('output/results_',runs{i}.name,'.mat');

    % If specified, import weights from weights file (tsv)
    if ~isempty(weightsPerReactionFile)
        weightsPerReactionFullPath = strcat(baseDirectory,'input/',weightsPerReactionFile);
        file_handle = fopen(weightsPerReactionFullPath);
        try
            u = textscan(file_handle,'%s\t%s');
        catch
            fprintf('File "%s" could not be found, exiting.\n',weightsPerReactionFullPath)
            return
        end
        weightsPerReaction.rxns = {};
        weightsPerReaction.weights = {};
        for k = 1:length(u{1})
            weightsPerReaction.rxns{k} = u{1}{k};
            weightsPerReaction.weights{k} = str2double(u{2}{k});
        end
        fclose(file_handle); 
    else
        weightsPerReaction = [];
    end
    

    cnt = 1;
%     cnt
%     run_idx
%     i
%     runs{i}
%     runs{i}.name
    Stats{cnt,run_idx} = runs{i}.name;cnt = cnt+1;
    % fastGapFill
    epsilon = 1e-4;
    tic; [AddedRxns] = fastGapFill(consistMatricesSUX,epsilon,weights,weightsPerReaction);
    tgap=toc;
    Stats{cnt,run_idx} = num2str(length(AddedRxns.rxns));cnt = cnt+1;
    save(SFilename);

    % Postprocessing
    [AddedRxnsExtended] = postProcessGapFillSolutions(AddedRxns,model,BlockedRxns,0);
    %clear AddedRxns;

    Stats{cnt,run_idx} = num2str(AddedRxnsExtended.Stats.metabolicSol);cnt = cnt+1;
    Stats{cnt,run_idx} = num2str(AddedRxnsExtended.Stats.transportSol);cnt = cnt+1;
    Stats{cnt,run_idx} = num2str(AddedRxnsExtended.Stats.exchangeSol);cnt = cnt+1;
    Stats{cnt,run_idx} = num2str(tgap);cnt = cnt+1;
    Stats{cnt,run_idx} = SFilename;cnt = cnt+1;
    clear a b

    % Reaction List
    col = 1;
    RxnList={};
    RxnList{1,col}=SFilename;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.rxns; col = col + 1;
    RxnList{1,col}=SFilename;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.rxnFormula; col = col + 1;
    RxnList{1,col}=SFilename;RxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.subSystem; col = col + 1;

    save(SFilename);
    
%     [m,n] = size(Stats);
%     fprintf('Stats are %ix%i (rows x columns)\n',m,n)
    
%     for j = 1:m
%         fprintf('%s: %f\n',Stats{j,1},Stats{j,run_idx});
%     end

end
    
%% Clean up

cd(oldDirectory)
clear run_directory model_file db_file dictionary_file;
