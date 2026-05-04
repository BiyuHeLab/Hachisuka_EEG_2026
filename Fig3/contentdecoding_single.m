%% EEG timeseries classification for content decoding (super/basic categories)
% This script not only balances the classes of interest, but
% the exemplars that make up the classes.

% Last updated Feb 6 2026, Ayaka Hachisuka

trialtype = 'ALL';
numCores = feature('numcores')-2;

conditionsLIST = {'Normal','Masked','LSF'};
superorbasicLIST = {'Super','Basic1','Basic2'};

run_pca = 0;
time_radius = 0;        
trial_bin_num = 4;
numrepeats = 300;

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');

for c = 1:length(conditionsLIST)
    condition = conditionsLIST{c};
for sb = 1:length(superorbasicLIST)

superorbasic = superorbasicLIST{sb};

if contains(superorbasic,'Super')
    categorylevel = 'Super';
elseif contains(superorbasic,'Basic')
    categorylevel = 'Basic';
elseif contains(superorbasic,'RecUnrec')
    categorylevel = 'RecUnrec';
end

%% get timelock data in CoSMoMVPA format
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');

% set configuration
config=cosmo_config();
config.data_path = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_input/']);

data_path=fullfile(config.data_path);
if contains(superorbasic,'Super')
    index2label={'Animate','Inanimate'}; % 1=pre, 2=peri/post
elseif contains(superorbasic,'Basic1')
    index2label={'Cat','Dog'};
elseif contains (superorbasic, 'Basic2')
    index2label={'Car','Truck'};
elseif contains (superorbasic, 'RecUnrec')
    index2label={'Rec','Unrec'};
end

for s = 1:33
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

%% Pre-allocate partitions: stratified CV

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
    % if testsize == 0; testsize = 1; end

    discardsub = 0;
    notenoughtrial = 0;

    if testsize == 0
        % trial_bin_num = 2;
        clearvars exempAall exempBall partitions
        exempAall = unique_values(ismember(unique_values,possiblexemp1));
        exempBall = unique_values(ismember(unique_values,possiblexemp2));
        
        if isempty(exempAall) || isempty(exempBall)
            disp('Not enough unique exemplars, discarding subject')
            discardsub = 1;
            break;
        else
        test_trials = []; test_trial_tempA = []; test_trial_tempB = [];
        sampling = 0;
        exempAall_orig = exempAall; exempBall_orig = exempBall;
        while sampling < (numtrials/2)*0.1 %total trials to be left out is 10% of all trials (divide by two, because a pair of exemplars are left out).

            % Reset exemp values when they are empty:
            if isempty(exempAall); exempAall = exempAall_orig; end
            if isempty(exempBall); exempBall = exempBall_orig; end

            % Randomly choose an exemplar from class A & B:

            % If there is only one exemplar, do not randomly sample:
            if length(exempAall) > 1; exempA = randsample(exempAall,1); else; exempA = exempAall; end
            if length(exempBall) > 1; exempB = randsample(exempBall,1); else; exempB = exempBall; end

            idxA = find(ds_tl2.sa.trialinfo == exempA);
            test_trial_tempA = randsample(idxA,1); %randomly sample 1 for class A

            idxB = find(ds_tl2.sa.trialinfo == exempB);
            test_trial_tempB = randsample(idxB,1); %randomly sample 1 for class B

            test_trials = [test_trials test_trial_tempA' test_trial_tempB'];

            % Remove the sampled exemplars:
            exempAall(exempAall == exempA) = [];
            exempBall(exempBall == exempB) = [];

            sampling = sampling + 1;
        end
        test_index = ismember(allindex,test_trials);
        train_index = ~test_index;
        partitions.train_indices{1} = allindex(train_index);
        partitions.test_indices{1} = allindex(test_index);
        end

    else
    % trial_bin_num = 4;
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

    end
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

if notenoughtrial == 0 && discardsub == 0
sl_map_samples = nan(numrepeats,241);
partitions = cell(numrepeats,1);
sl_map = cell(numrepeats,1);
parfor (rr = 1:numrepeats,numCores)
% for rr = 1:numrepeats
    ds_tl3 = ds_sel3{rr};
    partitions{rr}.test_indices{1} = partitionsall{rr}.test_indices{1};
    partitions{rr}.train_indices{1} = partitionsall{rr}.train_indices{1};
    bal_partitions=cosmo_balance_partitions(partitions{rr}, ds_tl3);
    
    % in addition give a label to each trial
    ds_tl3.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl3.sa.targets));
    
    % just to check everything is ok
    cosmo_check_dataset(ds_tl3);
    
    average_train_args_cell={{'average_train_count',trial_bin_num,...            %  { 2
                                    'average_train_resamplings',1}}; %  { 2
    n_average_train_args=numel(average_train_args_cell);
    
    chantypes = {'eeg'};
    
    nchantypes=length(chantypes);
    
    npartitions=numel(bal_partitions);
    fprintf('There are %d partitions\n', numel(bal_partitions.train_indices));
    fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                            bal_partitions.train_indices)));
    fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                            bal_partitions.test_indices)));
    measure=@cosmo_crossvalidation_measure; %@cosmo_crossvalidation_measure doesn't give you weights vector
    measure_args=struct();
    measure_args.classifier=@cosmo_classify_lda3;
    measure_args.partitions=bal_partitions;
    measure_args.output = 'balanced_accuracy';
    % measure_args.nproc = 38;
    
    if run_pca == 1
        measure_args.pca_explained_ratio = 0.95;
    end
    % add the options to average samples to the measure arguments.
    % (if no averaging is desired, this step can be left out.)
    if ~isempty(average_train_args_cell)
        average_train_args=average_train_args_cell{1};
        measure_args=cosmo_structjoin(measure_args, average_train_args);
    end
                          
    nbrhood=cosmo_interval_neighborhood(ds_tl3,'time',...
                                            'radius',time_radius);
    sl_map1{rr}=cosmo_searchlight(ds_tl3,nbrhood,measure,measure_args);
    sl_map_samples(rr,:) = sl_map1{rr}.samples; %save every repeat accuracy
    % end
end
    sl_map = sl_map{1};
    sl_map.samples = nanmean(sl_map_samples); %average across all repeats
    
    if ~isempty(sl_map)
            sl_map = sl_map1{1};
            samplesall = nan(numrepeats,241);
            weightsall = nan(120,241,numrepeats);
            weight_vec = nan(120,241);
            activation_pattern_all = nan(120,241,numrepeats);
            for r = 1:numrepeats
                for t = 1:length(sl_map1{r}.hyperplane)
                    temp = cell2mat(sl_map1{r}.hyperplane{t});
                    weight_vec(:,t) = temp'; %temp(2,:) - temp(1,:);
                end
                sl_map1{r}.hyperplane = [];
                sl_map1{r}.hyperplane = weight_vec;

                samplesall(r,:) = sl_map1{r}.samples;
                weightsall(:,:,r) = sl_map1{r}.hyperplane;
                for t = 1:length(sl_map1{r}.activation_pattern)
                    temp_pattern(:,t) = cell2mat(sl_map1{r}.activation_pattern{t});
                end
                activation_pattern_all(:,:,r) = temp_pattern;
            end
            sl_map.samples = nanmean(samplesall);
            sl_map.hyperplane = nanmean(weightsall,3);
            sl_map.activation_pattern = nanmean(activation_pattern_all,3);

    f = figure('visible','off'); set(gcf,'Color','w');
    hold on
    time_values=sl_map.a.fdim.values{1}; % first dim (channels got nuked)
    plot(time_values,sl_map.samples,'LineWidth',1.5,'Color','b');
    xline(0, '--', 'LineWidth', 1, 'Color', 'k');
    yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
    xlim([min(time_values),max(time_values)]);
    ylabel('Balanced Accuracy (chance=.5)');
    xlabel('time');

    descr=sprintf('Subject %s', num2str(s));
    title(descr);

    % save subject-level fig
    fig_save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechan',superorbasic,'_LDAba_',trialtype,'Trials/figs/']);
    if ~exist(fig_save_dir, 'dir'); mkdir(fig_save_dir); end
    filename = sprintf('sub%d.jpg',s);
    fullpath = fullfile(fig_save_dir, filename);
    saveas(gcf, fullpath);
    close all

    % Save data output files:
    save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechan',superorbasic,'_LDAba_',trialtype,'Trials/']);
    if ~exist(save_dir, 'dir'); mkdir(save_dir); end
    save([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput'],'sl_map','-v7.3');
    end
end
end
end
end
disp('Done!');
% end

