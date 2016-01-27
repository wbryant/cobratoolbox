function [AddedRxns] = submitFastGapFill(inputsFile,modelFileIn,dbFileIn,dictionaryFileIn,workspaceFileIn,weightsPerRxnFile,forceRerun,epsilon,blackList)
%% function [AddedRxns] = submitFastGapFill(modelFile,dbFile,dictionaryFile,workspaceFile,paramsFile)
%
% A convenience function for running both prepareFastGapFill and
% fastGapFill, allowing all files to be optionally specified and not 
% requiring prepareFastGapFill to be rerun (if relevant variables are saved
% in a specified workspace file) if it has already been run.
%
% This function performs two runs of fastGapFill, with weights as specified
% in the 'Define Runs' section of the code.  The 'runs' variable can be
% appended to or changed to enable the running of fastGapFill for different
% combinations of weights (and different weightsPerRxnFile files).
%
% INPUT (ALL optional, defaults to E. coli model iAF1260 and example files)
% inputsFile          .tsv file containing input file locations (can be used 
%                        instead of specifying input files in function call)
% modelFileIn         File containing model, either .mat or .xml
% dbFileIn            File containing universal database
% dictionaryFileIn    File containing db metabolite IDs and model
%                        counterparts, either .xls or .tsv
% workspaceFileIn     File for storing prepareFastGapFill results 
%                        (default: examples/defaultWorkspace.mat)
% weightsPerRxnFile   File containing individual weights for reactions 
%                        (default: empty)
% forceRerun          Rerun prepareFastGapFill even if it has already been run?
%                        (default: false)
% epsilon             fastCore parameter (default: 1e-4)
% blackList           List of excluded universal DB reactions 
%                        (default: none)
%
% OUTPUT
% AddedRxns           Reactions that have been added from UX matrix to S
%
% Jan 2016
% Will Bryant


%% Preparation - load data

% Suppress load warnings and deal with in the function
warning('off','MATLAB:load:variableNotFound')

% Default files
runFile = which('submitFastGapFill');
runDirCell = regexp(runFile,'(.+/)[^/]+$','tokens');
runDir = runDirCell{1}{1};

modelFile = strcat(runDir,'examples/iAF1260.mat');
dbFile = strcat(runDir,'AuxillaryFiles/reaction.lst');
dictionaryFile = strcat(runDir,'AuxillaryFiles/KEGG_dictionary.xls');
workspaceFile = strcat(runDir,'examples/defaultWorkspace.mat');

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

if ~exist('forceRerun','var') || isempty(forceRerun)
    forceRerun=false;
end

if ~exist('epsilon','var') || isempty(epsilon)
    epsilon=1e-4;
end

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
    fprintf('prepareGapFill already run, and forceRerun set to "false", prepareGapFill will not be rerun\n');
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
    
    % Save data
    save(workspaceFile,'consistModel','consistMatricesSUX','BlockedRxns','prepStats','model');
    fprintf('prepareGapFill finished\n')
end


prompt='Proceed with fastGapFill? [y]/n: ';
str = input(prompt,'s');
if ~isempty(str) && (str == 'n')
    return
end

fprintf('Running fastGapFill ...\n')
% define weights for reactions to be added - the lower the weight the
% higher the priority


%% Define Runs
% WARNING: weights should not be 1, as this impacts on the algorithm
% performance


if ~exist('weightsPerRxnFile','var') || isempty(weightsPerRxnFile)
    weightsPerRxnFile='';
end
runs = {};

% RUN 1 
run.weightsPerReactionFile = weightsPerRxnFile;
run.weights.MetabolicRxns = 0.1; % Kegg metabolic reactions
run.weights.ExchangeRxns = 0.5; % Exchange reactions
run.weights.TransportRxns = 10; % Transport reactions
run.name = 'initial';
runs = [runs; run];

% RUN 2
run.weightsPerReactionFile = weightsPerRxnFile;
run.weights.MetabolicRxns = 0.1; % Kegg metabolic reactions
run.weights.ExchangeRxns = 0.9; % Exchange reactions
run.weights.TransportRxns = 40; % Transport reactions
run.name = 'high_weight';
runs = [runs; run];

clear run

%% Run fastGapFill

cnt = 1;
clear Stats;
Stats{cnt,1} = 'Run name';cnt = cnt+1;
Stats{cnt,1} = 'Number of added reactions (all)';cnt = cnt+1;
Stats{cnt,1} = 'Number of added metabolic reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added transport reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added exchange reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Time fastGapFill';cnt = cnt+1;
Stats{cnt,1} = 'SFilename';cnt = cnt+1;

try
    a = load(workspaceFile,'consistModel','consistMatricesSUX','BlockedRxns','model');
    assert(length(fieldnames(a)) == 4);
    clear a;
    load(workspaceFile,'consistModel','consistMatricesSUX','BlockedRxns','model');
    fprintf('Variables loaded from .mat file\n');
    try
        b = load(workspaceFile,'prepStats');
        assert(~isempty(b))
        load(workspaceFile,'prepStats');
        clear b;
    catch
        fprintf('Warning: no stats for prepareFastGapfill\n');
    end
catch
    fprintf('Using preloaded variables for fastGapFill\n')
end

% [~,rxnSides] = regexp(rxnFormula,'(.+)<=+(.+)','match','tokens');
% wkspFolder = regexp(workspace,'(.+/)[^/]+$','tokens');



for i = 1:length(runs) 
    run_idx = i+1;
    fprintf('\nRun %i (%s)\n',i,runs{i}.name);
    weights = runs{i}.weights;
    weightsPerReactionFile = runs{i}.weightsPerReactionFile;
    SFilename = strcat(runDir,'examples/results_',runs{i}.name,'.mat');

    % If specified, import weights from weights file (tsv)
    if ~isempty(weightsPerReactionFile)
        file_handle = fopen(weightsPerReactionFile);
        try
            u = textscan(file_handle,'%s\t%s');
        catch
            fprintf('File "%s" could not be found, exiting\n',weightsPerReactionFile)
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
    Stats{cnt,run_idx} = runs{i}.name;cnt = cnt+1;
    
    % fastGapFill
    tic; 
    [AddedRxns] = fastGapFill(consistMatricesSUX,epsilon,weights,weightsPerReaction);
    tgap=toc;
    Stats{cnt,run_idx} = num2str(length(AddedRxns.rxns));cnt = cnt+1;
    save(SFilename);

    % Postprocessing
    [AddedRxnsExtended] = postProcessGapFillSolutions(AddedRxns,model,BlockedRxns,0);
    
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
end
    
%% Clean up
clear run_directory model_file db_file dictionary_file;
