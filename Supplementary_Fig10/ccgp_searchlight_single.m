%% CCGP Searchlight Subject-level analysis

% searchlight decoding analysis at the channel level.
% Doing this for cross-category decoding (training & testing on different
% basic categories) because the Haufe method relies on the covariance
% matrix of the training data, and here, the test results are highly
% critical!
% March 6 2025, Ayaka Hachisuka

%% Run searchlight
conditionLIST = {'Normal','Masked','LSF'};

for cc = 1:length(conditionLIST)
condition = conditionLIST{cc};
numCores = feature('numcores')-2;

addpath('/gpfs/data/helab/ah5385/EEG/Preprocessing');
addpath('/gpfs/data/helab/ah5385/EEG/fieldtrip-release/');

run_pca = 0;
trialtype = 'ALL';
time_radius = 0;        
trial_bin_num = 4;
numsubjects = 33;
numrepeats = 100;

addpath('/gpfs/data/helab/ah5385/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/gpfs/data/helab/ah5385/EEG/toolboxes/CoSMoMVPA/external/');

superorbasic = 'BasicABCD';
pairingList = {'ACBD','ADBC'};
analysistype = 'ALL'; %'rec' or 'unrec' or 'rec2' or 'unrec2'
analysis = sprintf(['supervsbasic',analysistype]);

% rec/unrec is the boundary between PAS0 vs. PAS1-3; rec2/unrec2 is the
% boundary between PAS0-1 vs. PAS2-3.

%% get timelock data in CoSMoMVPA format
addpath('/gpfs/data/helab/ah5385/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/gpfs/data/helab/ah5385/EEG/toolboxes/CoSMoMVPA/external/');

% set configuration
config=cosmo_config();
config.data_path = fullfile(['/gpfs/data/helab/ah5385/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_input/']);

data_path=fullfile(config.data_path);
index2label = {'1','2'};

% reset citation list
cosmo_check_external('-tic');

for s = 1:numsubjects

for pair = 1:2

if pair == 1
    A = 1;
    B = 3;
    C = 2;
    D = 4;
    % train_modality = [1 3]; test_modality = [2 4]; 
    pairing = 'ACBD';
elseif pair == 2
    A = 1;
    B = 4;
    C = 2;
    D = 3;
    pairing = 'ADBC';
    % train_modality = [1 4]; test_modality = [2 3];
end
    
subjNum = s;
if length(num2str(subjNum)) == 1 || ...
        (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
    subjNumStr = ['0' num2str(subjNum)];
else
    subjNumStr = num2str(subjNum);
end

% % load data

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
addpath('/gpfs/data/helab/ah5385/EEG/toolboxes/CoSMoMVPA/mvpa/');
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
average_train_args_cell={{'average_train_count',trial_bin_num,...
                                'average_train_resamplings',1}};
n_average_train_args=numel(average_train_args_cell);

chantypes = {'eeg'};
chan_count=length(chantypes);

index2label = {'dog','cat','car','truck'};
ds_tl.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl.sa.targets));
ds_tl.sa.exemptrialinfo = exemptrialinfo; %data_tl.trialinfo;

ds_sel3 = cell(n_modalities,n_modalities);
bal_partitionsall = cell(n_modalities,n_modalities);
    for train_modality=1:n_modalities
        for test_modality=1:n_modalities
            ds_sel3{train_modality,test_modality} = cell(numrepeats,1);
            bal_partitionsall{train_modality,test_modality} = cell(numrepeats,1);

            % select data in train and test modality
            msk=cosmo_match(ds_tl.sa.modality,[train_modality test_modality]);
            ds_sel=cosmo_slice(ds_tl,msk);
            
            for repeat = 1:numrepeats

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
        end
    end

    %% Prepare layout
easycaplayout = load('/gpfs/data/helab/ah5385/EEG/fromThomas/Data/layout_acticap128chan.mat').lay_acticap128chan;

common_labels = intersect(easycaplayout.label,ds_sel3{1,1}{1,1}.a.fdim.values{1,1});
[~, idx] = ismember(common_labels, easycaplayout.label);
new_layout = struct();

new_layout.pos = easycaplayout.pos(idx, :);
new_layout.width = easycaplayout.width(idx);
new_layout.height = easycaplayout.height(idx);
new_layout.label = easycaplayout.label(idx);
new_layout.outline = easycaplayout.outline;
new_layout.mask = easycaplayout.mask;
new_layout.cfg = easycaplayout.cfg;
new_layout.cfg.channel = new_layout.label;

cfg.layout = new_layout;
cfg.method = 'triangulation';
lay = ft_prepare_neighbours(cfg)'; %cosmo_meeg_chan_neighbors should do the same as this funtion          

clearvars ds_tl ds_sel ds_sel2 partitions
disp('Ready to start parallel process...');
for tr=1:n_modalities
for te=1:n_modalities
    if notenoughtrial(tr,te) == 0
    ds_sel_temp = ds_sel3{tr,te};
    bal_part_temp = bal_partitionsall{tr,te};
    parfor (rr = 1:numrepeats,numCores)
    % for rr = 1:numrepeats
    ds_sel4 = ds_sel_temp{rr};
    bal_partitions = bal_part_temp{rr};
  
    chan_nbrhood=cosmo_meeg_chan_neighborhood(ds_sel4, lay);
    time_nbrhood=cosmo_interval_neighborhood(ds_sel4,'time',...
                                                'radius',time_radius);
    
    % cross neighborhoods for chan-time searchlight
    nbrhood=cosmo_cross_neighborhood(ds_sel4,{chan_nbrhood,...
                                            time_nbrhood});
    
    % print some info
    nbrhood_nfeatures=cellfun(@numel,nbrhood.neighbors);
    fprintf('Features have on average %.1f +/- %.1f neighbors\n', ...
                mean(nbrhood_nfeatures), std(nbrhood_nfeatures));

    % just to check everything is ok
    cosmo_check_dataset(ds_sel4);
    npartitions=numel(bal_partitions);

    % only keep features with at least 10 neighbors
    % center_ids=find(nbrhood_nfeatures>10);
    measure=@cosmo_crossvalidation_measure;

    measure_args=struct();

    average_train_args_cell={{'average_train_count',trial_bin_num,...            %  { 2
                                    'average_train_resamplings',1}}; %  { 2
    n_average_train_args=numel(average_train_args_cell);
    
    if ~isempty(average_train_args_cell)
        average_train_args=average_train_args_cell{1};
        measure_args=cosmo_structjoin(measure_args, average_train_args);
    end

    measure_args.partitions=bal_partitions;
    measure_args.classifier=@cosmo_classify_lda3;
    measure_args.output = 'balanced_accuracy';
    measure_args.normalize = 'zscore';
    % run searchlight
    sl_tl_ds1{rr}=cosmo_searchlight(ds_sel4,nbrhood,measure,measure_args,'progress',[]);
    end
    if ~isempty(sl_tl_ds1)
        sl_tl_ds = sl_tl_ds1{1};
        for r = 1:numrepeats
            sl_tl_ds_samples(r,:) = sl_tl_ds1{r}.samples;
            % ds_searchlight_result_hyp(r,:) = ds_searchlight_result1{r}.hypall;
        end
        sl_tl_ds.samples = nanmean(sl_tl_ds_samples);
        
        savedir = fullfile(['/gpfs/data/helab/ah5385/EEG/CCGP_inputs/resultswith200Hzdata/Searchlight/',condition,'_BasicABCD', analysistype]);
        if ~exist(savedir,'dir'); mkdir(savedir); end
        save(fullfile(savedir, ['/sub', num2str(s), '_pair', num2str(pair), '_combo_', ...
            num2str(tr), num2str(te), '_',num2str(numrepeats),'repeats_sl_tl_ds']),'sl_tl_ds','-v7.3');
    end
    clearvars sl_tl_ds sl_tl_ds1

end
end
end
end
end
end
