%% MEEG timeseries classification

% A leave-one-trial out distance calculation.

% DECISION-to-BOUND analysis
% 36 exemplars total, 12 per condition.
% Final script for either SUPER or BASIC1/BASIC2.

% With 100 repeats

% Last updated April 23 2026, Ayaka Hachisuka

function distancetobound_contentdecoding_bytrial(trialtype)
conditionList = {'AllCond'};
num_subjects = 33;
subjectslist = setdiff(1:num_subjects,[22 32 33]);

for cc = 1:length(conditionList)
numrepeats = 100;
numCores = feature('numcores')-2;
run_pca = 0;
condition = conditionList{cc};

time_radius = 0;        
trial_bin_num = 1; %fixed value for this analysis
trial_resample = 1; %fixed

superorbasicList = {'Basic1','Basic2','Super'};

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing')

%% get timelock data in CoSMoMVPA format
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');

% set configuration
config=cosmo_config();
config.data_path = fullfile('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_inputv2/');
data_path=fullfile(config.data_path);

for ss = 1:length(superorbasicList)
    superorbasic = superorbasicList{ss};
if contains(superorbasic,'Super')
    index2label={'Animate','Inanimate'}; % 1=pre, 2=peri/post
    if strcmp(condition,'Normal')
        exemp_vals = 1:3:36;
    elseif strcmp(condition,'Masked')
        exemp_vals = 2:3:36;
    elseif strcmp(condition,'LSF')
        exemp_vals = 3:3:36;
    else
        exemp_vals = 1:36;
    end
elseif contains(superorbasic,'Basic1')
    index2label={'Cat','Dog'};
    if strcmp(condition,'Normal')
        exemp_vals = 1:3:18;
    elseif strcmp(condition,'Masked')
        exemp_vals = 2:3:18;
    elseif strcmp(condition,'LSF')
        exemp_vals = 3:3:18;
    else
        exemp_vals = 1:18;
    end
elseif contains (superorbasic, 'Basic2')
    index2label={'Car','Truck'};
    if strcmp(condition,'Normal')
        exemp_vals = 19:3:36;
    elseif strcmp(condition,'Masked')
        exemp_vals = 20:3:36;
    elseif strcmp(condition,'LSF')
        exemp_vals = 21:3:36;
    else
        exemp_vals = 19:36;
    end
elseif contains (superorbasic, 'RecUnrec')
    index2label={'Rec','Unrec'};
end
% reset citation list
cosmo_check_external('-tic');

for s = 1:length(subjectlist)

    subjNum = s;
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end

    data_fn=fullfile([data_path,'sub', subjNumStr, '_dataEEG_',trialtype,'.mat']);
    data_tl1=load(data_fn);
    data_tl = data_tl1.data3;
    trialnumbers = 1:length(data_tl.trialinfo);

    % convert to cosmomvpa struct
    addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
    ds_tl=cosmo_meeg_dataset(data_tl);
    
    %% Prepare MVPA
    ds_tl.sa.exemptrialinfo = data_tl.trialinfo;
    t = ds_tl.sa.exemptrialinfo;
    ds_tl.sa.targets = zeros(size(t,1),1);
    % set the target (trial condition)
    if strcmp(superorbasic,'Basic1') && strcmp(condition,'AllCond')
        ds_tl.sa.targets(t == 1 | t == 4 | t == 7) = 1; % dog (N)
        ds_tl.sa.targets(t == 10 | t == 13 | t == 16) = 2; % cat (N)
        ds_tl.sa.targets(t == 2 | t == 5 | t == 8) = 1; % dog (M)
        ds_tl.sa.targets(t == 11 | t == 14 | t == 17) = 2; % cat (M)
        ds_tl.sa.targets(t == 3 | t == 6 | t == 9) = 1; % dog (L)
        ds_tl.sa.targets(t == 12 | t == 15 | t == 18) = 2; % cat(L)
    elseif strcmp(superorbasic,'Basic1') && strcmp(condition,'Normal')
        ds_tl.sa.targets(t == 1 | t == 4 | t == 7) = 1; % dog (N)
        ds_tl.sa.targets(t == 10 | t == 13 | t == 16) = 2; % cat (N)
    elseif strcmp(superorbasic,'Basic1') && strcmp(condition,'Masked')
        ds_tl.sa.targets(t == 2 | t == 5 | t == 8) = 1; % dog (M)
        ds_tl.sa.targets(t == 11 | t == 14 | t == 17) = 2; % cat (M)
    elseif strcmp(superorbasic,'Basic1') && strcmp(condition,'LSF')
        ds_tl.sa.targets(t == 3 | t == 6 | t == 9) = 1; % dog (L)
        ds_tl.sa.targets(t == 12 | t == 15 | t == 18) = 2; % cat(L)
    elseif strcmp(superorbasic,'Basic2') && strcmp(condition,'AllCond')
        ds_tl.sa.targets(t == 19 | t == 22 | t == 25) = 1; % car (N)
        ds_tl.sa.targets(t == 28 | t == 31 | t == 34) = 2; % truck (N)
        ds_tl.sa.targets(t == 20 | t == 23 | t == 26) = 1; % car (M)
        ds_tl.sa.targets(t == 29 | t == 32 | t == 35) = 2; % truck (M)
        ds_tl.sa.targets(t == 21 | t == 24 | t == 27) = 1; % car (L)
        ds_tl.sa.targets(t == 30 | t == 33 | t == 36) = 2; % truck (L)
    elseif strcmp(superorbasic,'Basic2') && strcmp(condition,'Normal')
        ds_tl.sa.targets(t == 19 | t == 22 | t == 25) = 1; % car (N)
        ds_tl.sa.targets(t == 28 | t == 31 | t == 34) = 2; % truck (N)
    elseif strcmp(superorbasic,'Basic2') && strcmp(condition,'Masked')
        ds_tl.sa.targets(t == 20 | t == 23 | t == 26) = 1; % car (M)
        ds_tl.sa.targets(t == 29 | t == 32 | t == 35) = 2; % truck (M)
    elseif strcmp(superorbasic,'Basic2') && strcmp(condition,'LSF')
        ds_tl.sa.targets(t == 21 | t == 24 | t == 27) = 1; % car (L)
        ds_tl.sa.targets(t == 30 | t == 33 | t == 36) = 2; % truck (L)
    elseif strcmp(superorbasic,'Super') && strcmp(condition,'AllCond')
        ds_tl.sa.targets(t >= 1 & t <= 18) = 1; 
        ds_tl.sa.targets(t >= 19 & t <= 36) = 2; 
    elseif strcmp(superorbasic,'Super') && strcmp(condition,'Normal')
        idx = ismember(t, 1:3:18);
        idx2 = ismember(t, 19:3:36);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
    elseif strcmp(superorbasic,'Super') && strcmp(condition,'Masked')
        idx = ismember(t, 2:3:18);
        idx2 = ismember(t, 1:3:36);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
    elseif strcmp(superorbasic,'Super') && strcmp(condition,'LSF')
        idx = ismember(t, 1:3:36);
        idx2 = ismember(t, 1:3:36);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
    end

    ds_tl.sa.origtrialinfo = data_tl.origtrialinfo;
    nantrials = find(ds_tl.sa.targets == 0);
    ds_tl.sa.targets(nantrials) = [];
    ds_tl.sa.trialinfo(nantrials,:,:) = [];
    ds_tl.samples(nantrials,:) = [];
    ds_tl.sa.origtrialinfo(nantrials,:) = [];
    ds_tl.sa.exemptrialinfo(nantrials,:) = [];

    data_tl.trialinfo(nantrials,:) = [];
    trialnumbers(nantrials) = [];

    % in addition give a label to each trial
    ds_tl.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl.sa.targets));
    
    % just to check everything is ok
    cosmo_check_dataset(ds_tl);
    
    average_train_args_cell={{'average_train_count',trial_bin_num,...            %  { 2
                                    'average_train_resamplings',1}}; %  { 2
    n_average_train_args=numel(average_train_args_cell);
    
    chantypes = {'eeg'};
    clearvars partitions

    ds_tl.sa.chunks = (1:size(ds_tl.samples,1))';
    bal_partitions = cell(numrepeats,1);
    ds_sel2 = cell(numrepeats,1);
    whichtrial_all = cell(numrepeats,1);
    for rr = 1:numrepeats   
        foldcounter=0;
        %balance exemplar count:
        allexempvals = ds_tl.sa.exemptrialinfo;
        mask = ismember(allexempvals, exemp_vals); %exemp_vals is 1:36 for Super, 1:18 for Basic1, 19:36 for Basic2
        filtered_exemp = allexempvals(mask);
    
        [unique_vals, ~, indices] = unique(filtered_exemp);
        exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
        min_exemp_count = min(exemp_counts); %number of exemplar with lowest repeat, we are going to match counts for all other exemplars
        balanced_indices_msk = false(size(ds_tl.sa.exemptrialinfo)); %logical index mask
    
        % Loop over each unique exemplar and reduce trial count to minimum exemplar count:
        for i = unique_vals'
            exemp_ind = find(ds_tl.sa.exemptrialinfo == i);
            rand_selection = randsample(exemp_ind, min_exemp_count); %sample without replacement to minimize value 
            balanced_indices_msk(rand_selection) = true;
        end
    
        %remove all the other trials
        ds_sel = cosmo_slice(ds_tl,balanced_indices_msk);
        allindex = 1:size(ds_sel.samples,1); % label trials 1 to N.
        whichtrial_all{rr} = find(balanced_indices_msk); %vector of all TEST trial indices
        for fold = allindex
            foldcounter=foldcounter+1;
            test_index = fold;
            train_index = setdiff(allindex,test_index);

            balexemp_target1_ind = train_index(ds_sel.sa.targets(train_index) == 1);
            balexemp_target2_ind = train_index(ds_sel.sa.targets(train_index) == 2);
            min_target = min(numel(balexemp_target1_ind),numel(balexemp_target2_ind));
            balexemp_target1_ind = randsample(balexemp_target1_ind,min_target);
            balexemp_target2_ind = randsample(balexemp_target2_ind,min_target);
            
            bal_partitions{rr}.train_indices{fold} = [balexemp_target1_ind balexemp_target2_ind];
            bal_partitions{rr}.test_indices{fold} = test_index;
        end
            ds_sel2{rr} = ds_sel;
    end

    measure=@cosmo_crossvalidation_measure;
    nbrhood=cosmo_interval_neighborhood(ds_tl,'time',...
                                            'radius',time_radius);

    distforrepeat = cell(1,numrepeats); 
    samplesforrepeat = cell(1,numrepeats); 
    origtrialforrepeat = cell(1,numrepeats); 
    parfor (rr = 1:numrepeats, numCores)
    % for rr = 1:numrepeats

        measure_args=struct();
        measure_args.classifier=@cosmo_classify_lda4; %lda4 version is needed for this analysis!
        
        measure_args.output = 'fold_dist';

        if run_pca == 1
            measure_args.pca_explained_ratio = 0.95;
        end
    
        if ~isempty(average_train_args_cell)
            average_train_args=average_train_args_cell{1};
            measure_args=cosmo_structjoin(measure_args, average_train_args);
        end
        temp1 = nan(864,241);
        temp2 = nan(864,241);
        temp3 = nan(864,25);

        ds_sel3 = ds_sel2{rr};
        measure_args.partitions=bal_partitions{rr};
        sl_map=cosmo_searchlight(ds_sel3,nbrhood,measure,measure_args);
        whichtrial = whichtrial_all{rr}; %ds_sel3.trialnumbers;
        
        if iscell(whichtrial)
            whichtrial = cell2mat(whichtrial);
        end
        sl_map.trialnumbers = whichtrial;

        temp1(whichtrial,:) = sl_map.dist;
        temp2(whichtrial,:) = sl_map.samples;
        temp3(whichtrial,:) = ds_sel3.sa.origtrialinfo;
    
        distforrepeat{rr} = temp1;
        samplesforrepeat{rr}  = temp2;
        origtrialforrepeat{rr}  = temp3;
    end

    sl_map2.dist = [];
    sl_map2.samples = [];
    sl_map2.origtrialinfo = [];

    distfinal1 = cat(3, distforrepeat{:});
    sampfinal = cat(3, samplesforrepeat{:});
    trialinfofinal = cat(3, origtrialforrepeat{:});

    sl_map2.dist = nanmean(distfinal1,3);
    sl_map2.samples = nanmean(sampfinal,3);
    sl_map2.origtrialinfo = nanmean(trialinfofinal,3);

    % Save data output files:
    save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/DistancetoBound_distbytrial_eucdistv4_withrepeat/',superorbasic,'/',condition,'/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);

    if ~exist(save_dir, 'dir'); mkdir(save_dir); end
    save([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput'],'sl_map2','-v7.3');
    % save([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_ds_tl'],'ds_tl','-v7.3');
    clc;
    end
end
end
end
