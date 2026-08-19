%% MEEG time-by-time transfer classification

% function C2F_content_timegen_BP_final4(trialtype)

% numCores = feature('numcores')-2;
numCores = 20;
% conditionLIST = {'Normal','Masked','LSF'};
% superorbasicLIST = {'Super','Basic1','Basic2'}; %{'RecUnrec','RecUnrec_01vs23'};
conditionLIST = {'LSF'};
superorbasicLIST = {'Super'}; %{'RecUnrec','RecUnrec_01vs23'};
% run_pca = 0;
trialtype = 'ALL';
numrepeats = 300;
trial_bin_num = 4;

for c = 1:length(conditionLIST)
    condition = conditionLIST{c};
for sb = 1:length(superorbasicLIST)
superorbasic = superorbasicLIST{sb};

if contains(superorbasic,'Super')
    categorylevel = 'Super';
elseif contains(superorbasic,'Basic')
    categorylevel = 'Basic';
elseif contains(superorbasic,'RecUnrec')
    categorylevel = 'RecUnrec';
end

%only take the conditions of interest:
if contains(condition, 'Normal')
    cond1 = 1; cond2 = 4;
elseif contains(condition, 'Masked')
    cond1 = 2; cond2 = 5;
elseif contains(condition, 'LSF')
    cond1 = 3; cond2 = 6;
end

if contains(superorbasic,'Super')
    index2label={'Animate','Inanimate'}; % 1=pre, 2=peri/post
elseif contains(superorbasic,'Basic1')
    index2label={'Cat','Dog'};
elseif contains (superorbasic, 'Basic2')
    index2label={'Car','Truck'};
%     index2label={'Dog','Cat','Car','Truck'}; % 1=pre, 2=peri/post
elseif contains (superorbasic, 'RecUnrec')
    index2label={'Rec','Unrec'};
end

% set configuration
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');

config=cosmo_config();
config.data_path = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_input/']);

data_path=fullfile(config.data_path);

for s = 26:31 %1:31
    disp(s);
    subjNum = s;
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end

% % load data
data_fn=fullfile([data_path,'sub', subjNumStr, '_dataEEG_',trialtype,'.mat']);
%trialinfo is 1-36. Conditions are: Normal, Masked, LSF...iterative.
%Exemplars are cat1 cat1 cat1, cat2 cat2 cat2...truck3 truck3 truck3.

data_tl1=load(data_fn);
data_tl = data_tl1.data3;
clearvars data_tl1

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

if strcmp(superorbasic, 'Super')
indices = find(data_tl.trialinfo ~= cond1 & data_tl.trialinfo ~= cond2 & data_tl.trialinfo ~= cond3 & data_tl.trialinfo ~= cond4 ...
    & data_tl.trialinfo ~= cond5 & data_tl.trialinfo ~= cond6 & data_tl.trialinfo ~= cond7 & data_tl.trialinfo ~= cond8 ...
    & data_tl.trialinfo ~= cond9 & data_tl.trialinfo ~= cond10 & data_tl.trialinfo ~= cond11 & data_tl.trialinfo ~= cond12);
    possiblexemp1 = 1:6; possiblexemp2 = 7:12;
elseif strcmp(superorbasic,'Basic1')
    indices = find(data_tl.trialinfo ~= cond1 & data_tl.trialinfo ~= cond2 & data_tl.trialinfo ~= cond3 & data_tl.trialinfo ~= cond4 ...
    & data_tl.trialinfo ~= cond5 & data_tl.trialinfo ~= cond6);
    possiblexemp1 = 1:3; possiblexemp2 = 4:6;
elseif strcmp(superorbasic,'Basic2')
    indices = find(data_tl.trialinfo ~= cond7 & data_tl.trialinfo ~= cond8 ...
    & data_tl.trialinfo ~= cond9 & data_tl.trialinfo ~= cond10 & data_tl.trialinfo ~= cond11 & data_tl.trialinfo ~= cond12);
    possiblexemp1 = 1:3; possiblexemp2 = 4:6;
end

data_tl.trialinfo(indices) = [];
data_tl.trial(indices,:,:) = [];
data_tl.sampleinfo(indices,:) = [];

% FIND all the unique exempar numbers and re-number them, 1-12.
[unique_values, ~, new_indices] = unique(data_tl.trialinfo);
temptrialinfo = new_indices;
data_tl.trialinfo = temptrialinfo;

% convert to cosmomvpa struct
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
ds_tl=cosmo_meeg_dataset(data_tl);

%% Prepare MVPA
ds_tl.sa.chunks = (1:size(ds_tl.samples,1))';
t = data_tl.trialinfo;
clearvars data_tl
temptrialinfo = t;
if strcmp(superorbasic,'Super')
    ds_tl.sa.targets( t == 1 | t == 2 | t == 3 | t == 4 | t == 5 | t == 6) = 1; % animate
    ds_tl.sa.targets(t == 7 | t == 8 | t == 9 | t == 10 | t == 11 | t == 12) = 2; % inanimate
    ds_tl.sa.targets = ds_tl.sa.targets';
elseif contains(superorbasic,'Basic')
    ds_tl.sa.targets(t == 1 | t == 2 | t == 3) = 1; % Dog or Car
    ds_tl.sa.targets(t == 4 | t == 5 | t == 6) = 2; % Cat or Truck
    ds_tl.sa.targets = ds_tl.sa.targets';
end
% samplesize = length(ds_tl.a.fdim.values{2,1});

%% Pre-allocate partitions:

clearvars partitionsall data_tl
repeat = 0;
totalrepeat = numrepeats;
while repeat < totalrepeat
    repeat = repeat + 1;
    % disp(repeat);
    %balance exemplar count:
    [unique_values, ~, indices] = unique(ds_tl.sa.trialinfo);
    exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
    min_exemp_count = min(exemp_counts); %number of exemplar with lowest repeat, we are going to match counts for all other exemplars
    balanced_indices_msk = false(size(ds_tl.sa.trialinfo)); %logical index mask

    % Loop over each unique exemplar and reduce trial count to
    % minimum exemplar count:
    for i = 1:length(unique_values)
        exemp_ind = find(ds_tl.sa.trialinfo == unique_values(i));
        rand_selection = randsample(exemp_ind, min_exemp_count); %sample without replacement to minimize value 
        balanced_indices_msk(rand_selection) = true;
    end

    %remove all the other trials
    ds_tl2 = cosmo_slice(ds_tl,balanced_indices_msk);
    allindex = (1:size(ds_tl2.samples,1))'; %vector of all trial indices
    [unique_values, ~, indices] = unique(ds_tl2.sa.trialinfo); %unique exemplars present in dataset

    %Sample ~10% of allindex, with the restriction that each exemplar needs
    %to be sampled at least once.
    numtrials = sum(ds_tl2.sa.targets); %should be the same for 1 or 2
    numexemplars = length(unique_values);
    testsize = floor((numtrials*0.5)/numexemplars);
    if testsize == 0; testsize = 1; end
    clearvars partitions bal_partitions
    ds_tl2.sa.chunks = (1:size(ds_tl2.samples,1))';
    ds_tl2.sa.exemptrialinfo = ds_tl2.sa.trialinfo;
    test_trials = []; test_trial_temp = [];
    for exemp = unique_values'
        % Randomly extract two of each exemplar (~10% of trials).
        idx = find(ds_tl2.sa.exemptrialinfo == exemp);
        test_trial_temp = randsample(idx,testsize,false); %randomly sample 2 without replacement
        test_trials = [test_trials test_trial_temp'];
    end
    test_index = ismember(allindex,test_trials);
    train_index = ~test_index;
    partitions.train_indices{1} = allindex(train_index);
    partitions.test_indices{1} = allindex(test_index);
    ds_tl2.sa.chunks(partitions.train_indices{1}) = 1;
    ds_tl2.sa.chunks(partitions.test_indices{1}) = 2;
    % end
    %

    % if there are not enough training samples, discard the subject:
    num_target1_train = length(find(ds_tl2.sa.targets(allindex(partitions.train_indices{1})) == 1));
    num_target2_train = length(find(ds_tl2.sa.targets(allindex(partitions.train_indices{1})) == 2));
    %you need at least 2 samples per class in training:
    notenoughtrial = 0;

    if num_target1_train < ((2*trial_bin_num) + 1) || num_target2_train < ((2*trial_bin_num) + 1)
        disp('not enough trials! retrying...')
        repeat = repeat-1;
        totalrepeat = totalrepeat+1;
    else
        partitionsall{repeat}=partitions;
        ds_sel3{repeat} = ds_tl2;
    end
    if totalrepeat > 500
        disp('Not enough trials, discarding subject.');
        notenoughtrial = 1;
        break;
    end

end

% just to check everything is ok
cosmo_check_dataset(ds_tl);

% cdt_ds_samples = nan(14161,100);
cdt_ds1 = cell(numrepeats,1);
for rr = 1:numrepeats

ds_tl = ds_sel3{rr};

nchunks=2; % two chunks are required for this analysis
ds_tl.sa.chunks=cosmo_chunkize(ds_tl,nchunks); %found out this is unnecessary

ds_tl.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl.sa.targets));

%% time-by-time generalization on magneto- and gradio-meters seperately
% compute and plot accuracies for magnetometers and gradiometers separately
chan_types = {'eeg'};
nchan_types=numel(chan_types);
average_train_args_cell={{'average_train_count',trial_bin_num,...            %  { 2
                                    'average_train_resamplings',1}}; %  { 2
n_average_train_args=numel(average_train_args_cell);

% set arguments for the cosmo_dim_generalization_measure
measure_args=struct();
if ~isempty(average_train_args_cell)
    average_train_args=average_train_args_cell{1};
    measure_args=cosmo_structjoin(measure_args, average_train_args);
end
% the cosmo_dim_generalization_measure requires that another
% measure (here: the crossvalidation measure) is specified. The
% specified measure is applied for each combination of time points
measure_args.measure=@cosmo_crossvalidation_measure;

% When used ordinary, the cosmo_crossvalidation_measure itself
% requires two arguments:
% - classifier (here: LDA)
% - partitions
% However, because the cosmo_dim_generalization_measure defines
% the partitions itself, they are not set here.
measure_args.classifier=@cosmo_classify_lda3;
measure_args.output = 'balanced_accuracy';
measure_args.normalize = 'zscore';
% define the dimension over which generalization takes place
measure_args.dimension='time';
measure_args.nproc = numCores;
% define the radius for the time dimension. Here not just a single
% time-point is used, but also the time-point before it and the time-point
% after it.
measure_args.radius=1;

% make 'time' a sample dimension
% (this necessary for cosmo_dim_generalization_measure)
ds_time=cosmo_dim_transpose(ds_tl,'time',1);

cosmo_disp(ds_time);

% run transfer across time with the searchlight neighborhood
%cdt_ds=cosmo_cartesian_dim_transfer(ds_time,'time',measure,...
%                        'args',measure_args);
cdt_ds1{rr}=cosmo_dim_generalization_measure(ds_time,measure_args);
% cdt_ds_samples = cdt_ds1{rr}.samples;
end

if ~isempty(cdt_ds1)
    cdt_ds = cdt_ds1{1};
    for r = 1:numrepeats
        cdt_ds_samples(r,:) = cdt_ds1{r}.samples;
        % ds_searchlight_result_hyp(r,:) = ds_searchlight_result1{r}.hypall;
    end
    cdt_ds.samples = nanmean(cdt_ds_samples);
end

% : sa.labels has 57121 values in dimension 1, expected 1 (maybe the data was intended to be transposed?)
% cdt_ds.sa.labels = cdt_ds1{1}.sa.labels';
% unflatten the data to get train_time x test_time matrix
% [data, labels, values]=cosmo_unflatten(cdt_ds,1);
% dgm_ds_cell{k} = cdt_ds;

if ~isempty(cdt_ds1)
savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/Time_generalization_matrices/Timeradius_is_1/AllTrials/',superorbasic, '/', condition,'/']);
if ~exist(savedir); mkdir(savedir); end
% save([savedir,'/allsub_dgm_ds_cell'],'dgm_ds_cell');
save([savedir,'sub',num2str(subjNumStr),'_',condition,'_cdt_ds'],'cdt_ds');
end

end
end
end
% end

% We have now the results from cosmo_dim_generalization_measure in
% dgm_ds_cell. The next step is to change the dimensions so that
% the dimensions become feature dimensions. (this makes it possible
% to cluster the data.)
% group_cell=cell(n,1);
% for k=1:n
%     dgm_ds=dgm_ds_cell{k};
% 
%     % make train_time and test_time a feature dimension
%     ds=cosmo_dim_transpose(dgm_ds,{'train_time', 'test_time' },2);
% 
%     % for one-sample t-test
%     ds.sa.targets=1;
% 
%     % each participant is independent
%     ds.sa.chunks=k;
%     ds.a.fdim.values{1,1} = round(ds.a.fdim.values{1,1},3); %round time variable, because there is a miniscule jitter in the "0" value, creating problems with data stacking across subs.
%     ds.a.fdim.values{2,1} = round(ds.a.fdim.values{2,1},3);
%     group_cell{k}=ds;
% end
% 
% group_ds=cosmo_stack(group_cell);
% save([savedir,'/allsub_group_ds'],'group_ds');
% end
% end
