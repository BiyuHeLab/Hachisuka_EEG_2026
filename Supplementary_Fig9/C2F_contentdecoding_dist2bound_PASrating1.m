%% MEEG timeseries classification

% DECISION-to-BOUND analysis
% 36 exemplars total, 12 per condition.
% Final script for either SUPER or BASIC1/BASIC2.

% With 100 repeats
% Last updated Feb 6 2026, Ayaka Hachisuka (ahachisu@gmail.com)

conditionList = {'AllCond'};
superorbasicList = {'Super','Basic1','Basic2'};

for cc = 1:length(conditionList)
	numrepeats = 100;
	numCores = feature('numcores')-2;
	run_pca = 0;
	condition = conditionList{cc};

	trialtype = 'ALL'; % trialtype = 'SuperCorrect';
	time_radius = 0;        
	trial_bin_num = 1; %fixed value for this analysis
	trial_resample = 1; %fixed
	num_subjects = 31;
	
	% superorbasicList = {'Super'};

	addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing')
	% paths_EEG = EEG_SetPaths;
	% vars_EEG = EEG_SetVars(paths_EEG);
	% vars = vars_EEG;

	%% get timelock data in CoSMoMVPA format
	addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
	addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');

	% set configuration
	config.data_path = fullfile('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/PASrating1/');

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

	for s = 1:num_subjects
        fprintf('Running Subject %d...\n', s);
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
	    % trialnumbers = 1:length(data_tl.trialinfo);
        trialnumbers = 1:length(data_tl.trialinfo(:,2));
	    origtrialcount = trialnumbers;

	    % convert to cosmomvpa struct
	    addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
	    ds_tl=cosmo_meeg_dataset(data_tl);
	    
	    %% Prepare MVPA
	    ds_tl.sa.exemptrialinfo = data_tl.trialinfo(:,2);
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
		idx2 = ismember(t, 20:3:36);
		ds_tl.sa.targets(idx) = 1;
		ds_tl.sa.targets(idx2) = 2;
	    elseif strcmp(superorbasic,'Super') && strcmp(condition,'LSF')
		idx = ismember(t, 3:3:36);
		idx2 = ismember(t, 21:3:36);
		ds_tl.sa.targets(idx) = 1;
		ds_tl.sa.targets(idx2) = 2;
	    end

	    nantrials = find(ds_tl.sa.targets == 0);

	    trialnumbers(nantrials) = [];
	    ds_tl.sa.targets(nantrials) = [];
	    ds_tl.sa.trialinfo(nantrials,:) = [];
	    ds_tl.samples(nantrials,:) = [];
	    ds_tl.sa.exemptrialinfo(nantrials,:) = [];

	    % in addition give a label to each trial
	    ds_tl.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl.sa.targets));
	    
	    % just to check everything is ok
	    cosmo_check_dataset(ds_tl);
	    
	    % separate by position: all top is 1, all bottom is 2, all right is 3,
	    % all left is 4

	    if contains(condition, 'Normal')
		cond1 = 1; cond2 = 4; cond3 = 7; cond4 = 10;
	    elseif contains(condition, 'Masked')
		cond1 = 2; cond2 = 5; cond3 = 8; cond4 = 11;
	    elseif contains(condition, 'LSF')
		cond1 = 3; cond2 = 6; cond3 = 9; cond4 = 12;
	    elseif contains(condition, 'AllCond')
		cond1 = [1 2 3]; cond2 = [4 5 6]; cond3 = [7 8 9]; cond4 = [10 11 12];
	    end

	    fourcondvals = unique(ds_tl.sa.trialinfo(:,1));
	    if ~contains(condition, 'AllCond')
		    ds_tl.sa.trialinfo(ds_tl.sa.trialinfo(:,1) == fourcondvals(1),1) = 1;
		    ds_tl.sa.trialinfo(ds_tl.sa.trialinfo(:,1) == fourcondvals(2),1) = 2;
		    ds_tl.sa.trialinfo(ds_tl.sa.trialinfo(:,1) == fourcondvals(3),1) = 3;
		    ds_tl.sa.trialinfo(ds_tl.sa.trialinfo(:,1) == fourcondvals(4),1) = 4;
	    else
		    ds_tl.sa.trialinfo(ismember(ds_tl.sa.trialinfo(:,1),cond1),1) = 1;
		    ds_tl.sa.trialinfo(ismember(ds_tl.sa.trialinfo(:,1),cond2),1) = 2;
		    ds_tl.sa.trialinfo(ismember(ds_tl.sa.trialinfo(:,1),cond3),1) = 3;
		    ds_tl.sa.trialinfo(ismember(ds_tl.sa.trialinfo(:,1),cond4),1) = 4;
	    end

	    origtrialinfo = ds_tl.sa.trialinfo;

	    allsamples = nan(numel(origtrialcount),241,numrepeats);
	    alldist = nan(numel(origtrialcount),241,numrepeats);
	    trialinfofinal = nan(numel(origtrialcount),numrepeats);
	    posinfofinal = nan(numel(origtrialcount),numrepeats);
	    

		% remove indices that are not the position of interest for this
		% loop
		%otherposvals = setdiff(1:4,posvalue1);
		%indices = find(ds_tl.sa.trialinfo(:,1) ~= otherposvals(1) & ds_tl.sa.trialinfo(:,1) ~= otherposvals(2) & ds_tl.sa.trialinfo(:,1) ~= otherposvals(3));
		%ind4pos = find(ds_tl.sa.trialinfo(:,1) ~= posvalue1);

		ds_tl2 = ds_tl;
		trialnumbers2 = trialnumbers;

		% just to check everything is ok (again)
		cosmo_check_dataset(ds_tl2);

		chantypes = {'eeg'};
	    
		bal_partitions = cell(numrepeats,1);
		clearvars partitions bal_partitions
		ds_tl2.sa.chunks = (1:size(ds_tl2.samples,1))';
		allidx = 1:size(ds_tl2.samples,1);
        trialnumbers_final = cell(numrepeats,1);
		for rr = 1:numrepeats
            foldcounter = 0;
		    for fold = 1:size(ds_tl2.samples,1)
		        clearvars balanced_indices_msk
		        train_idx = setdiff(allidx,fold);
		        exemps = ds_tl2.sa.exemptrialinfo(train_idx);
		        %balance exemplar count:
		        [unique_values, ~, indices] = unique(exemps);
		        % disp('original unique value length');
		        % disp(length(unique_values));
		        exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
		        min_exemp_count = min(exemp_counts); %number of exemplar with lowest repeat, we are going to match counts for all other exemplars
		        balanced_indices_msk = false(length(train_idx),1); %logical index mask
		    
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
                    val1_temp = unique_values(ismember(unique_values,[2:3:18]));
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
    
                    if size(ds_tl2.samples,1) > 1
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
                        partitions.train_indices{foldcounter} = train_idx(balanced_indices_msk);
                        partitions.test_indices{foldcounter} = fold;
                    else
                        disp('not enough trials! discarding fold...');
                    end
                end

            if exist('partitions','var')
                bal_partitionsall{rr}=partitions; %cosmo_balance_partitions(partitions, ds_sel2);
                % in addition give a label to each trial
                ds_tl2.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl2.sa.targets));
                ds_sel{rr} = ds_tl2;
                trialnumbers_final{rr} = trialnumbers2;
            end
    		end
        end

        if exist('bal_partitionsall','var')
            if ~isempty(trialnumbers_final{rr})
		        measure=@cosmo_crossvalidation_measure;
		        nbrhood=cosmo_interval_neighborhood(ds_sel{rr},'time',...
		                                                'radius',time_radius);
	            
		        measure_args=struct();
		        measure_args.classifier=@cosmo_classify_lda4;
		        %measure_args.normalization = 'zscore';
		        measure_args.output = 'fold_dist';
	            
		        average_train_args_cell={{'average_train_count',trial_bin_num,...
		                        'average_train_resamplings',1}};
        
		        if ~isempty(average_train_args_cell)
		            average_train_args=average_train_args_cell{1};
		            measure_args=cosmo_structjoin(measure_args, average_train_args);
		        end
	            
		        for rr = 1:numrepeats
		            % if there are not enough training samples, discard the subject:
		            num_target1_train = length(find(ds_sel{rr}.sa.targets(bal_partitionsall{rr}.train_indices{1}) == 1));
		            num_target2_train = length(find(ds_sel{rr}.sa.targets(bal_partitionsall{rr}.train_indices{1}) == 2));
        
		            if num_target1_train < ((2*trial_bin_num) + 1) || num_target2_train < ((2*trial_bin_num) + 1)
		                disp('not enough trials! discarding subject...');
		            else
		                whichtrial = trialnumbers_final{rr};
		                measure_args.partitions=bal_partitionsall{rr};
                        % sl_map=cosmo_searchlight_dist(ds_sel{rr}, nbrhood,measure,measure_args);
		                sl_map=cosmo_searchlight_dist(ds_sel{rr}, nbrhood,measure,measure_args,'nproc',numCores, 'progress',[]);
		                allsamples(whichtrial,:,rr) = sl_map.samples;
		                alldist(whichtrial,:,rr) = sl_map.dist;
		                trialinfofinal(whichtrial,rr) = ds_sel{rr}.sa.exemptrialinfo;
                    end
                end
            end
        end

        sl_map2.dist = [];
        sl_map2.samples = [];
        sl_map2.origtrialinfo = [];
    
        sl_map2.dist = nanmean(alldist,3);
        sl_map2.samples = nanmean(allsamples,3);
        sl_map2.origtrialinfo = nanmean(trialinfofinal,3);
    
        % Save data output files: do this after collecting distance data for each
        % of the 4 positions
        % Save data output files:
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/',superorbasic,'/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
        if ~exist(save_dir, 'dir'); mkdir(save_dir); end
        save([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput'],'sl_map2','-v7.3');
        fprintf('Saving Subject #%d! \n',s);
end
end
