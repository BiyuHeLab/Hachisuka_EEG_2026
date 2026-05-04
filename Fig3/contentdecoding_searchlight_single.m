%% MEEG time-lock searchlight
% searchlight decoding analysis at the channel level.
% Last updated Feb 6 2026, Ayaka Hachisuka

function contentdecoding_searchlight_batch(condition)

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/');

superorbasicLIST = {'Super','Basic1','Basic2'};
numCores = feature('numcores')-2;

run_pca = 0;

% condition = 'Normal';
% superorbasic = 'Super';
trialtype = 'ALL';
time_radius = 0;        
trial_bin_num = 4;
numsubjects = 33;
numrepeats = 100;

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');

%% get timelock data in CoSMoMVPA format

% set configuration
config=cosmo_config();
config.data_path = fullfile('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_input/');
data_path=fullfile(config.data_path);

for sb = 1:length(superorbasicLIST)
    superorbasic = superorbasicLIST{sb};

if contains(superorbasic,'Super')
    index2label={'Animate','Inanimate'}; % 1=pre, 2=peri/post
elseif contains(superorbasic,'Basic1')
    index2label={'Cat','Dog'};
elseif contains (superorbasic, 'Basic2')
    index2label={'Car','Truck'};
elseif contains(superorbasic,'RecUnrec')
    index2label={'Recognized','Unrecognized'};
end

for s = 1:numsubjects
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
temptrialinfo = t;
if strcmp(superorbasic,'Super')
    ds_tl.sa.targets( t == 1 | t == 2 | t == 3 | t == 4 | t == 5 | t == 6) = 1; % animate
    ds_tl.sa.targets(t == 7 | t == 8 | t == 9 | t == 10 | t == 11 | t == 12) = 2; % inanimate
    ds_tl.sa.targets = ds_tl.sa.targets';
    possiblexemp1 = 1:6;
    possiblexemp2 = 7:12;
elseif contains(superorbasic,'Basic')
    ds_tl.sa.targets(t == 1 | t == 2 | t == 3) = 1; % Dog or Car
    ds_tl.sa.targets(t == 4 | t == 5 | t == 6) = 2; % Cat or Truck
    ds_tl.sa.targets = ds_tl.sa.targets';
    possiblexemp1 = 1:3;
    possiblexemp2 = 4:6;
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
    allindex = 1:size(ds_tl2.samples,1)'; %vector of all trial indices
    [unique_values, ~, indices] = unique(ds_tl2.sa.trialinfo); %unique exemplars present in dataset

    %Sample ~10% of allindex, with the restriction that each exemplar needs
    %to be sampled at least once.
    numtrials = sum(ds_tl2.sa.targets); %should be the same for 1 or 2
    numexemplars = length(unique_values);
    testsize = floor((numtrials*0.1)/numexemplars);
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

clearvars ds_tl ds_tl2
clearvars partitions sl_map_samples sl_map
%%
% define neighborhood parameters for each dimension
chan_type = 'eeg';

easycaplayout = load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/layout_acticap128chan.mat').lay_acticap128chan;

common_labels = intersect(easycaplayout.label,ds_sel3{1}.a.fdim.values{1});
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

if notenoughtrial == 0
partitions = cell(numrepeats,1);
sl_tl_ds1 = cell(numrepeats,1);
% parfor (rr = 1:numrepeats,numCores)
for rr = 1:numrepeats
    ds_tl3 = ds_sel3{rr};
    partitions{rr}.test_indices{1} = partitionsall{rr}.test_indices{1};
    partitions{rr}.train_indices{1} = partitionsall{rr}.train_indices{1};
    bal_partitions=cosmo_balance_partitions(partitions{rr}, ds_tl3);
      
    chan_nbrhood=cosmo_meeg_chan_neighborhood(ds_tl3, lay);
    time_nbrhood=cosmo_interval_neighborhood(ds_tl3,'time',...
                                                'radius',time_radius);
    
    % cross neighborhoods for chan-time searchlight
    nbrhood=cosmo_cross_neighborhood(ds_tl3,{chan_nbrhood,...
                                            time_nbrhood});
    
    % print some info
    nbrhood_nfeatures=cellfun(@numel,nbrhood.neighbors);
    fprintf('Features have on average %.1f +/- %.1f neighbors\n', ...
                mean(nbrhood_nfeatures), std(nbrhood_nfeatures));

    % in addition give a label to each trial
    ds_tl3.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl3.sa.targets));
    
    % just to check everything is ok
    cosmo_check_dataset(ds_tl3);
    npartitions=numel(bal_partitions);
    fprintf('There are %d partitions\n', numel(bal_partitions.train_indices));
    fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                            bal_partitions.train_indices)));
    fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                            bal_partitions.test_indices)));

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

    % run searchlight
    sl_tl_ds1{rr}=cosmo_searchlight(ds_tl3,nbrhood,measure,measure_args,'progress',[]);
end
end
if ~isempty(sl_tl_ds1)
sl_tl_ds = sl_tl_ds1{1};
for r = 1:numrepeats
    sl_tl_ds_samples(r,:) = sl_tl_ds1{r}.samples;
    % ds_searchlight_result_hyp(r,:) = ds_searchlight_result1{r}.hypall;
end
sl_tl_ds.samples = nanmean(sl_tl_ds_samples);

save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/Searchlight/noPCA_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan',superorbasic,'_LDA_AllTrials/']);
if ~exist(save_dir, 'dir'); mkdir(save_dir); end
save([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'sl_tl_ds'],'sl_tl_ds','-v7.3');
end
end
end
end
