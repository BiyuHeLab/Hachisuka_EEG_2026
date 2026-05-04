%% CCGP Analysis: Training and testing across the SUPERORDINATE Boundary.
% Decoding analysis pipeline.
% For debugging, change numrepeats to 3 and batchsize to 1. Batch size is
% to run analysis in "batches" to reduce computational load.

% Final script was run on a high-performing cluster with 40 cores and 150GB
% memory; takes approximately 18 hours.

% Train on Dog vs. Car, Test on Cat vs. Truck.
% Train on Dog vs. Truck, Test on Cat vs. Car.
% Train on Cat vs. Truck, Test on Dog vs. Car.
% Train on Cat vs. Car, Test on Dog vs. Truck.

%Last updated June 20 2025, Ayaka Hachisuka

function ccgp_supervsbasic_batchv2(analysistype)
num_subjects = 33;
time_radius = 0;
numrepeats = 100;
batchsize = 10;
trial_bin_num = 4;
conditionLIST = {'Normal','Masked','LSF'};
numCores = feature('numcores')-2;

for ii = 1:length(conditionLIST)
condition = conditionLIST{ii};
superorbasic = 'BasicABCD';
pairingList = {'ACBD','ADBC'};
analysis = sprintf(['supervsbasic',analysistype]);

%% get timelock data in CoSMoMVPA format
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');

% set configuration
config=cosmo_config();
config.data_path = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_input/']);

data_path=fullfile(config.data_path);
index2label = {'1','2'};

% reset citation list
cosmo_check_external('-tic');

for s = 1:num_subjects

for pair = 1:2

if pair == 1
    A = 1;
    B = 3;
    C = 2;
    D = 4;
    pairing = 'ACBD';
elseif pair == 2
    A = 1;
    B = 4;
    C = 2;
    D = 3;
    pairing = 'ADBC';
end
    
subjNum = s;
if length(num2str(subjNum)) == 1 || ...
        (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
    subjNumStr = ['0' num2str(subjNum)];
else
    subjNumStr = num2str(subjNum);
end

%load data

data_fn=fullfile([data_path,'sub', subjNumStr, '_dataEEG_',analysistype,'.mat']);

data_tl1=load(data_fn);
data_tl = data_tl1.data3;

%only take the conditions of interest:
if contains(condition, 'Normal')
    cond1 = 1; cond2 = 4; cond3 = 7; cond4 = 10;
    cond5 = 13; cond6 = 16; cond7 = 19; cond8 = 22;
    cond9 = 25; cond10 = 28; cond11 = 31; cond12 = 34;
elseif contains(condition, 'Masked')
    cond1 = 2; cond2 = 5; cond3 = 8; cond4 = 11;
    cond5 = 14; cond6 = 17; cond7 = 20; cond8 = 23;
    cond9 = 26; cond10 = 29; cond11 = 32; cond12 = 35;
elseif contains(condition, 'LSF')
    cond1 = 3; cond2 = 6; cond3 = 9; cond4 = 12;
    cond5 = 15; cond6 = 18; cond7 = 21; cond8 = 24;
    cond9 = 27; cond10 = 30; cond11 = 33; cond12 = 36;
end

indices = find(data_tl.trialinfo ~= cond1 & data_tl.trialinfo ~= cond2 & data_tl.trialinfo ~= cond3 & data_tl.trialinfo ~= cond4 ...
    & data_tl.trialinfo ~= cond5 & data_tl.trialinfo ~= cond6 & data_tl.trialinfo ~= cond7 & data_tl.trialinfo ~= cond8 ...
    & data_tl.trialinfo ~= cond9 & data_tl.trialinfo ~= cond10 & data_tl.trialinfo ~= cond11 & data_tl.trialinfo ~= cond12);

data_tl.trialinfo(indices) = [];
data_tl.trial(indices,:,:) = [];
data_tl.sampleinfo(indices,:) = [];

% FIND all the unique exempar numbers and re-number them, 1-12.
[unique_values, ~, new_indices] = unique(data_tl.trialinfo);
temptrialinfo = new_indices;
data_tl.trialinfo = temptrialinfo;

exemptrialinfo = data_tl.trialinfo; % for trialinfo assignment, below.

% convert to cosmomvpa struct
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
ds_tl=cosmo_meeg_dataset(data_tl);
clearvars data_tl

ds_tl.sa.trialinfo(exemptrialinfo == 1 | exemptrialinfo == 2 | exemptrialinfo == 3) = 1; % dog, animate
ds_tl.sa.trialinfo(exemptrialinfo == 4 | exemptrialinfo == 5 | exemptrialinfo == 6) = 2; % cat, animate
ds_tl.sa.trialinfo(exemptrialinfo == 7 | exemptrialinfo == 8 | exemptrialinfo == 9) = 3; % car, inanimate
ds_tl.sa.trialinfo(exemptrialinfo == 10 | exemptrialinfo == 11 | exemptrialinfo == 12) = 4; % truck, inanimate

for t = 1:size(ds_tl.samples,1)
    if ds_tl.sa.trialinfo(t) == A 
        ds_tl.sa.targets(t,:) = 1;
    elseif ds_tl.sa.trialinfo(t) == B
        ds_tl.sa.targets(t,:) = 2;
    elseif ds_tl.sa.trialinfo(t) == C
        ds_tl.sa.targets(t,:) = 1;
    elseif ds_tl.sa.trialinfo(t) == D
        ds_tl.sa.targets(t,:) = 2;
    end
end

n_modalities=2;

for t = 1:size(ds_tl.samples,1)
    if ds_tl.sa.trialinfo(t) == A || ds_tl.sa.trialinfo(t) == B
        ds_tl.sa.modality(t,:) = 1;
    elseif ds_tl.sa.trialinfo(t) == C || ds_tl.sa.trialinfo(t) == D
        ds_tl.sa.modality(t,:) = 2;
    end
end  

ds_tl.sa.chunks = (1:size(ds_tl.samples,1))';

%% Prepare MVPA

% just to check everything is ok
cosmo_check_dataset(ds_tl);

average_train_args_cell={{'average_train_count',trial_bin_num,...            %  { 2
                                'average_train_resamplings',1}}; %  { 2
n_average_train_args=numel(average_train_args_cell);

chantypes = {'eeg'};
chan_count=length(chantypes);

index2label = {'dog','cat','car','truck'};
ds_tl.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl.sa.targets));
ds_tl.sa.exemptrialinfo = exemptrialinfo; %data_tl.trialinfo;
%% do all combinations for training and test modalities

% do all combinations for training and test modalities
ds_sel3 = cell(n_modalities,n_modalities);
bal_partitionsall = cell(n_modalities,n_modalities);
    for train_modality=1:n_modalities
        for test_modality=1:n_modalities
            ds_sel3{train_modality,test_modality} = cell(numrepeats,1);
            bal_partitionsall{train_modality,test_modality} = cell(numrepeats,1);

            for batch = 1:(numrepeats/batchsize)
            % Start and end indices for this batch
            startIdx = (batch - 1) * batchsize + 1;
            endIdx = batch * batchsize;
           
            % select data in train and test modality
            msk=cosmo_match(ds_tl.sa.modality,[train_modality test_modality]);
            ds_sel=cosmo_slice(ds_tl,msk);

            for repeat = startIdx:endIdx

            %balance exemplar count:
            [unique_values, ~, indices] = unique(ds_sel.sa.exemptrialinfo);
            exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
            min_exemp_count = min(exemp_counts); %number of exemplar with lowest repeat, we are going to match counts for all other exemplars
            balanced_indices_msk = false(size(ds_sel.sa.exemptrialinfo)); %logical index mask
        
            % Loop over each unique exemplar and reduce trial count to
            % minimum exemplar count:

            for i = unique_values'
                exemp_ind = find(ds_sel.sa.exemptrialinfo == i);
                rand_selection = randsample(exemp_ind, min_exemp_count); %sample without replacement to minimize value 
                balanced_indices_msk(rand_selection) = true;
            end
        
            %remove all the other trials
            ds_sel2 = cosmo_slice(ds_sel,balanced_indices_msk);
            allindex = 1:size(ds_sel2.samples,1); %vector of all trial indices
            [unique_values, ~, indices] = unique(ds_sel2.sa.exemptrialinfo); %unique exemplars present in dataset

            if train_modality~=test_modality
                clearvars partitions bal_partitions
                % Partition the balanced data:
                ds_sel2.sa.chunks = ds_sel2.sa.modality;
                partitions = cosmo_nchoosek_partitioner(ds_sel2,1,'modality',test_modality);
                % bal_partitionsall{train_modality,test_modality}{repeat}=cosmo_balance_partitions(partitions, ds_sel2);
                % ds_sel3{train_modality,test_modality}{repeat} = ds_sel2;
            else
                clearvars partitions bal_partitions
                ds_sel2.sa.chunks = (1:size(ds_sel2.samples,1))';
                allindex = 1:size(ds_sel2.samples,1);
                test_trials = []; test_trial_temp = [];
                for exemp = unique_values'
                    % Randomly extract two of each exemplar (~10% of trials).
                    idx = find(ds_sel2.sa.exemptrialinfo == exemp);
                    test_trial_temp = randsample(idx,1,false); %randomly sample 2 without replacement
                    test_trials = [test_trials test_trial_temp'];
                end
                test_index = ismember(allindex,test_trials);
                train_index = ~test_index;
                partitions.train_indices{1} = allindex(train_index);
                partitions.test_indices{1} = allindex(test_index);
            end
                % if there are not enough training samples, discard the subject:
                num_target1_train = length(find(ds_sel2.sa.targets(allindex(partitions.train_indices{1})) == 1));
                num_target2_train = length(find(ds_sel2.sa.targets(allindex(partitions.train_indices{1})) == 2));
                %you need at least 2 samples per class in training:
                notenoughtrial(train_modality,test_modality) = 0;

                if num_target1_train < ((2*trial_bin_num) + 1) || num_target2_train < ((2*trial_bin_num) + 1)
                    disp('Not enough trials...');
                    notenoughtrial(train_modality,test_modality) = 1;
                    break;
                else
                    bal_partitionsall{train_modality,test_modality}{repeat}=cosmo_balance_partitions(partitions, ds_sel2);
                    ds_sel3{train_modality,test_modality}{repeat} = ds_sel2;
                end
            end

    clearvars ds_sel ds_sel2 partitions
    disp('Starting parallel process...');
            if notenoughtrial(train_modality,test_modality) == 0
            ds_sel_temp = ds_sel3{train_modality,test_modality};
            bal_part_temp = bal_partitionsall{train_modality,test_modality};
            ds_searchlight_result_samples = nan(numrepeats,241);
            ds_searchlight_result1 = cell(numrepeats,1);
            % parfor (rr = 1:numrepeats,numCores)
            for rr = 1:numrepeats
            disp(startIdx);
            disp(endIdx);
            parfor (rr = startIdx:endIdx,numCores)
                ds_sel4 = ds_sel_temp{rr};
                bal_partitions = bal_part_temp{rr};
                opt=struct();
                opt.partitions=bal_partitions;
                opt.classifier=@cosmo_classify_lda3;
                opt.max_feature_count = 30000;
                opt.output = 'balanced_accuracy';

                if ~isempty(average_train_args_cell)
                    average_train_args=average_train_args_cell{1};
                    opt=cosmo_structjoin(opt, average_train_args);
                end

                measure=@cosmo_crossvalidation_measure;
    
                % Run the measure.
                nh=cosmo_interval_neighborhood(ds_sel4,'time','radius',time_radius);
                ds_searchlight_result1{rr}=cosmo_searchlight(ds_sel4,nh,measure,opt,'progress',[]); %searchlight option
                tempsamples(rr,:) = ds_searchlight_result1{rr}.samples; %save every repeat accuracy
                temphyperplane(rr,:) = ds_searchlight_result1{rr}.hyperplane;
                tempactpattern(rr,:) = ds_searchlight_result1{rr}.activation_pattern;
            % end
            end
            ds_samples{train_modality,test_modality}(startIdx:endIdx,:) = tempsamples(startIdx:endIdx,:);
            ds_hyperplane{train_modality,test_modality}(startIdx:endIdx,:) = temphyperplane(startIdx:endIdx,:);
            ds_actpattern{train_modality,test_modality}(startIdx:endIdx,:) = tempactpattern(startIdx:endIdx,:);
            end
            end
        end
    end
    for train_modality=1:n_modalities
    for test_modality=1:n_modalities
        if ~isempty(ds_searchlight_result1)
        ds_searchlight_result = ds_searchlight_result1{end};
        for t = 1:size(ds_hyperplane{train_modality,test_modality},2)
            temp_weight = ds_hyperplane{train_modality,test_modality}(:,t);
            temp_weight2 = cell2mat(cellfun(@(x) x{1}(:)', temp_weight, 'UniformOutput', false));
            weight_vec(:,t) = nanmean(temp_weight2)';

            temp_pattern = ds_actpattern{train_modality,test_modality}(:,t);
            temp_pattern2 = cell2mat(cellfun(@(x) x{1}(:)', temp_pattern, 'UniformOutput', false));
            act_pattern(:,t) = nanmean(temp_pattern2)';

            temp_samples = ds_samples{train_modality,test_modality}(:,t);
            samplesall(:,t) = nanmean(temp_samples)';
        end

        ds_searchlight_result.samples = samplesall;
        ds_searchlight_result.hyperplane = weight_vec;
        ds_searchlight_result.activation_pattern = act_pattern;

        savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/resultswith200Hzdata/',condition,'_BasicABCD', analysistype]);
        if ~exist(savedir,'dir'); mkdir(savedir); end
        fprintf('Saving data for pair %s combo %s %s\n', num2str(pair), num2str(train_modality), num2str(test_modality));
        save(fullfile(savedir, ['/sub', num2str(s), '_pair', num2str(pair), '_combo_', ...
            num2str(train_modality), num2str(test_modality), '_',num2str(numrepeats),'repeats_ds_searchlight_result']),'ds_searchlight_result','-v7.3');
        clearvars ds_searchlight_result
        end
        end
    end
end
end
end
end

