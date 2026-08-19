%% Position decoding
%normal, top - 1
%masked, top - 2
%lsf, top - 3
%normal, bottom - 4
%masked, bottom - 5
%lsf, bottom - 6
%normal, left - 7
%masked, left - 8
%lsf, left - 9 
%normal, right - 10
%masked, right - 11
%lsf, right - 12

clear
close all

conditionlist = {'Normal','Masked','LSF'}; %{'AllCond'};

superorbasicLIST = {'Super','Basic1','Basic2};
trialtype = 'ALL';

numrepeats = 300;
numCores = feature('numcores');

for ss = 1:length(superorbasicLIST)
superorbasic = superorbasicLIST{ss};

for cc = 1:length(conditionlist)

condition = conditionlist{cc};

run_pca = 0;
time_radius = 0;        
trial_bin_num = 4;
trial_resample = 1;
nchunks = 10;
nminval = 1;
num_subjects = 31;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing')
% paths_EEG = EEG_SetPaths;
% vars_EEG = EEG_SetVars(paths_EEG);
% vars = vars_EEG;

%% get timelock data in CoSMoMVPA format
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');

config=cosmo_config();
config.data_path = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/PositionDecoding/final_decoding_input/']);

data_path=fullfile(config.data_path);
index2label = {'Position1','Position2'};

% reset citation list
cosmo_check_external('-tic');

for s = 1:num_subjects
    fprintf('Running Subject %d...\n', s);
    ds_sel2 = cell(numrepeats,1);
    subjNum = s;
    fprintf('Running Subject %d...\n', s);
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end

    data_fn=fullfile([data_path,'sub', subjNumStr, '_dataEEG_',trialtype,'.mat']);
    data_tl1=load(data_fn);
    data_tl = data_tl1.data3;

    %only take the conditions of interest:
    if contains(condition, 'Normal')
        cond1 = 1; cond2 = 4; cond3 = 7; cond4 = 10;
    elseif contains(condition, 'Masked')
        cond1 = 2; cond2 = 5; cond3 = 8; cond4 = 11;
    elseif contains(condition, 'LSF')
        cond1 = 3; cond2 = 6; cond3 = 9; cond4 = 12;
    elseif contains(condition, 'AllCond')
        cond1 = [1 2 3]; cond2 = [4 5 6]; cond3 = [7 8 9]; cond4 = [10 11 12];
    end

    indices = find(~ismember(data_tl.trialinfo(:,1), cond1) & ~ismember(data_tl.trialinfo(:,1), cond2) & ~ismember(data_tl.trialinfo(:,1), cond3) & ~ismember(data_tl.trialinfo(:,1), cond4));
    data_tl.trialinfo(indices,:) = [];
    data_tl.trial(indices,:,:) = [];
    data_tl.sampleinfo(indices,:) = [];

    fourcondvals = unique(data_tl.trialinfo(:,1));
    if ~contains(condition, 'AllCond')
        data_tl.trialinfo(data_tl.trialinfo(:,1) == fourcondvals(1),1) = 1;
        data_tl.trialinfo(data_tl.trialinfo(:,1) == fourcondvals(2),1) = 2;
        data_tl.trialinfo(data_tl.trialinfo(:,1) == fourcondvals(3),1) = 3;
        data_tl.trialinfo(data_tl.trialinfo(:,1) == fourcondvals(4),1) = 4;
    else
        data_tl.trialinfo(ismember(data_tl.trialinfo(:,1),cond1),1) = 1;
        data_tl.trialinfo(ismember(data_tl.trialinfo(:,1),cond2),1) = 2;
        data_tl.trialinfo(ismember(data_tl.trialinfo(:,1),cond3),1) = 3;
        data_tl.trialinfo(ismember(data_tl.trialinfo(:,1),cond4),1) = 4;
    end
    % do the same for superordinate category distinction:
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
    
  if ~contains(condition, 'AllCond')
        if strcmp(superorbasic, 'Super')
        indices = find(data_tl.trialinfo(:,2) ~= cond1 & data_tl.trialinfo(:,2) ~= cond2 & data_tl.trialinfo(:,2)~= cond3 & data_tl.trialinfo(:,2)~= cond4 ...
            & data_tl.trialinfo(:,2)~= cond5 & data_tl.trialinfo(:,2)~= cond6 & data_tl.trialinfo(:,2)~= cond7 & data_tl.trialinfo(:,2)~= cond8 ...
            & data_tl.trialinfo(:,2)~= cond9 & data_tl.trialinfo(:,2)~= cond10 & data_tl.trialinfo(:,2)~= cond11 & data_tl.trialinfo(:,2)~= cond12);
            exemp_vals = 1:12;
        elseif strcmp(superorbasic,'Basic1')
            indices = find(data_tl.trialinfo(:,2)~= cond1 & data_tl.trialinfo(:,2)~= cond2 & data_tl.trialinfo(:,2)~= cond3 & data_tl.trialinfo(:,2)~= cond4 ...
            & data_tl.trialinfo(:,2)~= cond5 & data_tl.trialinfo(:,2)~= cond6);
            exemp_vals = 1:6;
        elseif strcmp(superorbasic,'Basic2')
            indices = find(data_tl.trialinfo(:,2)~= cond7 & data_tl.trialinfo(:,2)~= cond8 ...
            & data_tl.trialinfo(:,2)~= cond9 & data_tl.trialinfo(:,2)~= cond10 & data_tl.trialinfo(:,2)~= cond11 & data_tl.trialinfo(:,2)~= cond12);
            exemp_vals = 7:12;
        end
        
        data_tl.trialinfo(indices,:) = [];
        data_tl.trial(indices,:,:) = [];
        data_tl.sampleinfo(indices,:) = [];
    else
      % exemp_vals is 1:36 for Super, 1:18 for Basic1, 19:36 for Basic2
        if strcmp(superorbasic, 'Super')
            exemp_vals = 1:36;
        elseif strcmp(superorbasic,'Basic1')
            exemp_vals = 1:18;
        elseif strcmp(superorbasic,'Basic2')
            exemp_vals = 19:36;
        end
    end
    
    % FIND all the unique exempar numbers and re-number them, 1-12.
    % [unique_values, ~, new_indices] = unique(data_tl.trialinfo(:,2));

    % new_indices(indices) = [];
    % data_tl.trial(indices,:,:) = [];
    % data_tl.sampleinfo(indices,:) = [];

    % convert to cosmomvpa struct
    ds_tl=cosmo_meeg_dataset(data_tl);
    % ds_tl.sa.trialinfo(:,2) = new_indices;

     %% Prepare MVPA
    ds_tl.sa.exemptrialinfo = ds_tl.sa.trialinfo(:,2);
    ds_tl.sa.posvalues = ds_tl.sa.trialinfo(:,1);

    t = ds_tl.sa.exemptrialinfo;
    ds_tl.sa.targets = zeros(size(t,1),1);
    % set the target (trial condition)
    if strcmp(superorbasic,'Basic1') && strcmp(condition,'AllCond')
        % ds_tl.sa.targets(t == 1 | t == 4 | t == 7) = 1; % dog (N)
        % ds_tl.sa.targets(t == 10 | t == 13 | t == 16) = 2; % cat (N)
        % ds_tl.sa.targets(t == 2 | t == 5 | t == 8) = 1; % dog (M)
        % ds_tl.sa.targets(t == 11 | t == 14 | t == 17) = 2; % cat (M)
        % ds_tl.sa.targets(t == 3 | t == 6 | t == 9) = 1; % dog (L)
        % ds_tl.sa.targets(t == 12 | t == 15 | t == 18) = 2; % cat(L)
        idx = ismember(t, 1:9);
        idx2 = ismember(t, 10:18);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
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
        % ds_tl.sa.targets(t == 19 | t == 22 | t == 25) = 1; % car (N)
        % ds_tl.sa.targets(t == 28 | t == 31 | t == 34) = 2; % truck (N)
        % ds_tl.sa.targets(t == 20 | t == 23 | t == 26) = 1; % car (M)
        % ds_tl.sa.targets(t == 29 | t == 32 | t == 35) = 2; % truck (M)
        % ds_tl.sa.targets(t == 21 | t == 24 | t == 27) = 1; % car (L)
        % ds_tl.sa.targets(t == 30 | t == 33 | t == 36) = 2; % truck (L)
        idx = ismember(t, 19:27);
        idx2 = ismember(t, 28:36);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
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
        idx = ismember(t, 1:18);
        idx2 = ismember(t, 19:36);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
    elseif strcmp(superorbasic,'Super') && strcmp(condition,'Normal')
        idx = ismember(t, 1:3:18);
        idx2 = ismember(t, 19:3:36);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
    elseif strcmp(superorbasic,'Super') && strcmp(condition,'Masked')
        idx = ismember(t, 2:3:18);
        idx2 = ismember(t, 20:3:36);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
    elseif strcmp(superorbasic,'Super') && strcmp(condition,'LSF')
        idx = ismember(t, 3:3:18);
        idx2 = ismember(t, 21:3:36);
        ds_tl.sa.targets(idx) = 1;
        ds_tl.sa.targets(idx2) = 2;
    end

    nantrials = find(ds_tl.sa.targets == 0);

    ds_tl.sa.targets(nantrials) = [];
    ds_tl.sa.posvalues(nantrials,:) = [];
    ds_tl.sa.trialinfo(nantrials,:) = [];
    ds_tl.samples(nantrials,:) = [];
    ds_tl.sa.exemptrialinfo(nantrials,:) = [];

    % Split train/test by position: train on 3 positions, test on the 4th
    for posvalue1 = 1:4
        otherposvals = setdiff(1:4,posvalue1);
        train_indices = find(ds_tl.sa.trialinfo(:,1) ~= posvalue1);
        test_indices = find(ds_tl.sa.trialinfo(:,1) ~= otherposvals(1) & ds_tl.sa.trialinfo(:,1) ~= otherposvals(2) & ds_tl.sa.trialinfo(:,1) ~= otherposvals(3));

        ds_tl.sa.modality = nan(size(ds_tl.samples,1),1);
        ds_tl.sa.modality(train_indices) = 1;
        ds_tl.sa.modality(test_indices) = 2;

        ds_tl.sa.chunks = (1:size(ds_tl.samples,1))';

        t = ds_tl.sa.trialinfo(:,2);
        temptrialinfo = t;

        clearvars bal_partitions
        %remove all the other trials
        for rr = 1:numrepeats
            clearvars partitions
            foldcounter = 0;
            for fold = 1:length(test_indices)
                test_index = test_indices(fold);
                exemps = temptrialinfo(train_indices);
                
                %balance exemplar count:
                [unique_values, ~, indices] = unique(exemps);
                % disp('original unique value length');
                % disp(length(unique_values));
                exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
                min_exemp_count = min(exemp_counts); %number of exemplar with lowest repeat, we are going to match counts for all other exemplars
                balanced_indices_msk = false(length(train_indices),1); %logical index mask
            
                % make sure there is an equal number of animate/inanimate or
                % cat/dog, car/truck trials in the unique_values (became a
                % concern for instances of fewer trials like PAS = 0)
                % Count numbers from 1 to 25
                if contains(superorbasic,'Super') && contains(condition,'AllCond')
                    count1 = sum(ismember(unique_values,[1:18]));
                    count2 = sum(ismember(unique_values,[19:36]));
                    val1_temp = unique_values(ismember(unique_values,[1:18]));
                    val2_temp = unique_values(ismember(unique_values,[19:36]));
                elseif strcmp(superorbasic,'Basic1') && contains(condition,'AllCond')
                    count1 = sum(ismember(unique_values,[1:9]));
                    count2 = sum(ismember(unique_values,[10:18]));
                    val1_temp = unique_values(ismember(unique_values,[1:9]));
                    val2_temp = unique_values(ismember(unique_values,[10:18]));
                elseif strcmp(superorbasic,'Basic2') && contains(condition,'AllCond')
                    count1 = sum(ismember(unique_values,[19:27]));
                    count2 = sum(ismember(unique_values,[28:36]));
                    val1_temp = unique_values(ismember(unique_values,[19:27]));
                    val2_temp = unique_values(ismember(unique_values,[28:36]));
                elseif contains(superorbasic,'Super') && contains(condition,'Normal')
                    count1 = sum(ismember(unique_values,[1:3:18]));
                    count2 = sum(ismember(unique_values,[19:3:36]));
                    val1_temp = unique_values(ismember(unique_values,[1:3:18]));
                    val2_temp = unique_values(ismember(unique_values,[19:3:36]));
                elseif strcmp(superorbasic,'Basic1') && contains(condition,'Normal')
                    count1 = sum(ismember(unique_values,[1:3:9]));
                    count2 = sum(ismember(unique_values,[10:3:18]));
                    val1_temp = unique_values(ismember(unique_values,[1:3:9]));
                    val2_temp = unique_values(ismember(unique_values,[10:3:18]));
                elseif strcmp(superorbasic,'Basic2') && contains(condition,'Normal')
                    count1 = sum(ismember(unique_values,[19:3:27]));
                    count2 = sum(ismember(unique_values,[28:3:36]));
                    val1_temp = unique_values(ismember(unique_values,[19:3:27]));
                    val2_temp = unique_values(ismember(unique_values,[28:3:36]));
                elseif contains(superorbasic,'Super') && contains(condition,'Masked')
                    count1 = sum(ismember(unique_values,[2:3:18]));
                    count2 = sum(ismember(unique_values,[20:3:36]));
                    val1_temp = unique_values(ismember(unique_values,[2:3:18]));
                    val2_temp = unique_values(ismember(unique_values,[20:3:36]));
                elseif strcmp(superorbasic,'Basic1') && contains(condition,'Masked')
                    count1 = sum(ismember(unique_values,[2:3:9]));
                    count2 = sum(ismember(unique_values,[11:3:18]));
                    val1_temp = unique_values(ismember(unique_values,[2:3:9]));
                    val2_temp = unique_values(ismember(unique_values,[11:3:18]));
                elseif strcmp(superorbasic,'Basic2') && contains(condition,'Masked')
                    count1 = sum(ismember(unique_values,[20:3:27]));
                    count2 = sum(ismember(unique_values,[29:3:36]));
                    val1_temp = unique_values(ismember(unique_values,[20:3:27]));
                    val2_temp = unique_values(ismember(unique_values,[29:3:36]));
                elseif contains(superorbasic,'Super') && contains(condition,'LSF')
                    count1 = sum(ismember(unique_values,[3:3:18]));
                    count2 = sum(ismember(unique_values,[21:3:36]));
                    val1_temp = unique_values(ismember(unique_values,[3:3:18]));
                    val2_temp = unique_values(ismember(unique_values,[21:3:36]));
                elseif strcmp(superorbasic,'Basic1') && contains(condition,'LSF')
                    count1 = sum(ismember(unique_values,[3:3:9]));
                    count2 = sum(ismember(unique_values,[12:3:18]));
                    val1_temp = unique_values(ismember(unique_values,[3:3:9]));
                    val2_temp = unique_values(ismember(unique_values,[12:3:18]));
                elseif strcmp(superorbasic,'Basic2') && contains(condition,'LSF')
                    count1 = sum(ismember(unique_values,[21:3:27]));
                    count2 = sum(ismember(unique_values,[30:3:36]));
                    val1_temp = unique_values(ismember(unique_values,[21:3:27]));
                    val2_temp = unique_values(ismember(unique_values,[30:3:36]));
                end

                if isempty(val1_temp) || isempty(val2_temp)
                    continue;
                else
                    if isscalar(val1_temp)
                        val1 = val1_temp;
                    else
                        val1 = randsample(val1_temp,min(count1,count2));
                    end
                    if isscalar(val2_temp)
                        val2 = val2_temp;
                    else
                        val2 = randsample(val2_temp,min(count1,count2));
                    end
                    unique_values2 = [val1; val2];
                
                    % Loop over each unique exemplar and reduce trial count to
                    % minimum exemplar count:
                    % disp('new unique value length');
                    % disp(length(unique_values));
                    exemps2 = exemps(ismember(exemps, unique_values2));
                    [unique_values, ~, indices] = unique(exemps2);
                    exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
                    min_exemp_count = min(exemp_counts); %number of exemplar with lowest repeat, we are going to match counts for all other exemplars
    
                    if size(ds_tl.samples,1) > 1
                        for i = unique_values2'
                            exemp_ind = find(exemps == i);
                            if isscalar(exemp_ind)
                                rand_selection = exemp_ind;
                            else
                                rand_selection = randsample(exemp_ind, min_exemp_count); %sample without replacement to minimize value 
                            end
                            balanced_indices_msk(rand_selection) = true;
                        end
                        foldcounter = foldcounter + 1;
                        partitions.train_indices{foldcounter} = train_indices(balanced_indices_msk);
                        partitions.test_indices{foldcounter} = test_index;
                    else
                        disp('not enough trials! discarding fold...');
                    end
                end
            end

            if exist('partitions','var')
                bal_partitionsall{rr}=partitions; %cosmo_balance_partitions(partitions, ds_sel2);
                % in addition give a label to each trial
                ds_tl.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl.sa.targets));
                ds_sel2{rr} = ds_tl;
            end
        end
    
        if exist('bal_partitionsall','var')
            measure_args=struct();
            measure_args.classifier=@cosmo_classify_lda3;
            measure_args.output = 'balanced_accuracy';
            %measure_args.normalization = 'zscore';
            average_train_args_cell={{'average_train_count',trial_bin_num,...            %  { 2
                                            'average_train_resamplings',1}}; %  { 2
            n_average_train_args=numel(average_train_args_cell);
            
            if ~isempty(average_train_args_cell)
                average_train_args=average_train_args_cell{1};
                measure_args=cosmo_structjoin(measure_args, average_train_args);
            end
            chantypes = {'eeg'};
    
            for rr = 1:numrepeats
    
                measure=@cosmo_crossvalidation_measure;
                nbrhood=cosmo_interval_neighborhood(ds_sel2{rr},'time',...
                                                    'radius',time_radius);
    
                % if there are not enough training samples, discard the subject:
                num_target1_train = length(find(ds_sel2{rr}.sa.targets(bal_partitionsall{rr}.train_indices{1}) == 1));
                num_target2_train = length(find(ds_sel2{rr}.sa.targets(bal_partitionsall{rr}.train_indices{1}) == 2));
                if num_target1_train <= ((2*trial_bin_num) + 1) || num_target2_train <= ((2*trial_bin_num) + 1)
                    disp('not enough trials! discarding subject...');
                    notenoughtrials = 1;
                    allsamples(rr,:) = nan(1,241);
                else
                    measure_args.partitions=bal_partitionsall{rr};
                    sl_map=cosmo_searchlight(ds_sel2{rr},nbrhood,measure,measure_args,'nproc',numCores,'progress',[]);
                    allsamples(rr,:) = sl_map.samples;
                end
            end
    
	        sl_map.samples = nanmean(allsamples);
	        % Save data output files:
	        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Final_Resubmission_results_April16/Fig3/CategoryDecodingCrossPos_leaveonetrialout_300repeats/',condition,'/',superorbasic,'/35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), '_position',num2str(posvalue1),'_LDAba_AllTrials/']);
	        if ~exist(save_dir, 'dir'); mkdir(save_dir); end
	        save([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput'],'sl_map','-v7.3');
	        fprintf('Savaing data for subject %d...\n', s);
        end
    end
end
end
end


