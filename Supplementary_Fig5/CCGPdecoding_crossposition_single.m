%% CCGP Analysis: Training and testing across the SUPERORDINATE Boundary.
% Decoding analysis pipeline.

% Train on Dog vs. Car, Test on Cat vs. Truck.
% Train on Dog vs. Truck, Test on Cat vs. Car.
% Train on Cat vs. Truck, Test on Dog vs. Car.
% Train on Cat vs. Car, Test on Dog vs. Truck.

analysistype = 'ALL';
num_subjects = 31;
time_radius = 0;
numrepeats = 100;
batchsize = 10;
trial_bin_num = 1;
conditionLIST = {'Normal','Masked','LSF'};
numCores = feature('numcores');

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
config.data_path = fullfile('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/PositionDecoding/final_decoding_input/');

data_path=fullfile(config.data_path);
index2label = {'1','2'};

% reset citation list
cosmo_check_external('-tic');

for s = [18 19 20 21 22 23 24 29 30 31]
fprintf('Running Subject %d...\n', s);
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

%load data

data_fn=fullfile([data_path,'sub', subjNumStr, '_dataEEG_',analysistype,'.mat']);

data_tl1=load(data_fn);
data_tl = data_tl1.data3;

%only take the conditions of interest:
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

    % indices = find(data_tl.trialinfo(:,1) ~= cond1 & data_tl.trialinfo(:,1) ~= cond2 & data_tl.trialinfo(:,1) ~= cond3 & data_tl.trialinfo(:,1) ~= cond4);
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


% FIND all the unique exempar numbers and re-number them, 1-12.
[unique_values, ~, new_indices] = unique(data_tl.trialinfo(:,2));
data_tl.trialinfo(:,2) = new_indices;

exemptrialinfo = data_tl.trialinfo(:,2); % for trialinfo assignment, below.

% convert to cosmomvpa struct
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
ds_tl=cosmo_meeg_dataset(data_tl);
clearvars data_tl

ds_tl.sa.trialinfo(ismember(exemptrialinfo,1:3),2) = 1; % dog, animate
ds_tl.sa.trialinfo(ismember(exemptrialinfo,4:6),2) = 2; % cat, animate
ds_tl.sa.trialinfo(ismember(exemptrialinfo,7:9),2) = 3; % car, inanimate
ds_tl.sa.trialinfo(ismember(exemptrialinfo,10:12),2) = 4; % truck, inanimate

for t = 1:size(ds_tl.samples,1)
    if ds_tl.sa.trialinfo(t,2) == A 
        ds_tl.sa.targets(t,:) = 1;
    elseif ds_tl.sa.trialinfo(t,2) == B
        ds_tl.sa.targets(t,:) = 2;
    elseif ds_tl.sa.trialinfo(t,2) == C
        ds_tl.sa.targets(t,:) = 1;
    elseif ds_tl.sa.trialinfo(t,2) == D
        ds_tl.sa.targets(t,:) = 2;
    end
end

n_modalities=2;

for t = 1:size(ds_tl.samples,1)
    if ds_tl.sa.trialinfo(t,2) == A || ds_tl.sa.trialinfo(t,2) == B
        ds_tl.sa.modality(t,:) = 1;
    elseif ds_tl.sa.trialinfo(t,2) == C || ds_tl.sa.trialinfo(t,2) == D
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

            for posvalue1 = 1:4
                otherposvals = setdiff(1:4,posvalue1);
		        train_pos = find(ds_tl.sa.trialinfo(:,1) ~= posvalue1);
		        test_pos = find(ds_tl.sa.trialinfo(:,1) ~= otherposvals(1) & ds_tl.sa.trialinfo(:,1) ~= otherposvals(2) & ds_tl.sa.trialinfo(:,1) ~= otherposvals(3));

                %allindices = 1:size(ds_tl.sa.trialinfo,1);
                %position_mask = find(ds_tl.sa.trialinfo(:,1) == posvalue1);
                %indices = find(ds_tl.sa.trialinfo(:,1) ~= otherposvals(1) & ds_tl.sa.trialinfo(:,1) ~= otherposvals(2) & ds_tl.sa.trialinfo(:,1) ~= otherposvals(3));       
                %ds_sel2 = cosmo_slice(ds_tl,position_mask);
                ds_sel2 = ds_tl;
                
                for rr = 1:numrepeats
                clearvars partitions bal_partitions
                if train_modality~=test_modality
                    % Partition the balanced data:
                    % ds_sel2.sa.chunks = ds_sel2.sa.modality;
                    % partitions = cosmo_nchoosek_partitioner(ds_sel2,1,'modality',test_modality);
                    % bal_partitionsall{train_modality,test_modality}{repeat}=cosmo_balance_partitions(partitions, ds_sel2);
                    % ds_sel3{train_modality,test_modality}{repeat} = ds_sel2;
                    ds_sel2.sa.chunks = (1:size(ds_sel2.samples,1))';
                    allidx = 1:size(ds_sel2.samples,1);
                    foldcounter = 0;
                    for fold = 1:size(ds_sel2.samples,1)
                        clearvars balanced_indices_msk
                        train_idx = intersect(train_pos,find(ds_sel2.sa.modality == train_modality));
                        test_idx = intersect(test_pos,find(ds_sel2.sa.modality == test_modality));
                        exemps = ds_sel2.sa.exemptrialinfo(train_idx);
                    
                        %balance exemplar count:
                        [unique_values, ~, indices] = unique(exemps);
                        exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
                        min_exemp_count = min(exemp_counts); %number of exemplar with lowest repeat, we are going to match counts for all other exemplars
                        balanced_indices_msk = false(1,length(train_idx)); %logical index mask
                    
                        % make sure there is an equal number of animate/inanimate or
                        % cat/dog, car/truck trials in the unique_values (became a
                        % concern for instances of fewer trials like PAS = 0)
                        % Count numbers from 1 to 25
                        count1 = sum(unique_values >= 1 & unique_values <= 6);
                        count2 = sum(unique_values >= 7 & unique_values <= 12);
                        val1_temp = unique_values(unique_values <= 6);
                        val2_temp = unique_values(unique_values >= 7);
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
                            unique_values = [val1; val2];
                            % Loop over each unique exemplar and reduce trial count to
                            % minimum exemplar count:
                            for i = unique_values'
                                exemp_ind = find(exemps == i);
                                if size(exemp_ind,1) == 1
                                    rand_selection = exemp_ind;
                                else
                                    rand_selection = randsample(exemp_ind, min_exemp_count); %sample without replacement to minimize value 
                                end
                                balanced_indices_msk(rand_selection) = true;
                            end
                            foldcounter = foldcounter + 1;
                            partitions.train_indices{foldcounter} = train_idx(balanced_indices_msk);
                            partitions.test_indices{foldcounter} = test_idx;
                        end
                    end
                else
                    clearvars partitions bal_partitions
                    ds_sel2.sa.chunks = (1:size(ds_sel2.samples,1))';
                    allidx = 1:size(ds_sel2.samples,1);
                    foldcounter = 0;
                    for fold = test_pos
                        
                        clearvars balanced_indices_msk
                        % train_idx = setdiff(allidx,fold);
                        train_idx = train_pos;
                        exemps = ds_sel2.sa.exemptrialinfo(train_idx);
                        %balance exemplar count:
                        [unique_values, ~, indices] = unique(exemps);
                        exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
                        min_exemp_count = min(exemp_counts); %number of exemplar with lowest repeat, we are going to match counts for all other exemplars
                        balanced_indices_msk = false(length(train_idx),1); %logical index mask
                    
                        % make sure there is an equal number of animate/inanimate or
                        % cat/dog, car/truck trials in the unique_values (became a
                        % concern for instances of fewer trials like PAS = 0)
                        % Count numbers from 1 to 25
                        count1 = sum(unique_values >= 1 & unique_values <= 6);
                        count2 = sum(unique_values >= 7 & unique_values <= 12);
                        val1_temp = unique_values(unique_values <= 6);
                        val2_temp = unique_values(unique_values >= 7);
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
                            unique_values = [val1; val2];
                            % Loop over each unique exemplar and reduce trial count to
                            % minimum exemplar count:
                            for i = unique_values'
                                exemp_ind = find(exemps == i);
                                if size(exemp_ind,1) == 1
                                    rand_selection = exemp_ind;
                                else
                                    rand_selection = randsample(exemp_ind, min_exemp_count); %sample without replacement to minimize value 
                                end
                                balanced_indices_msk(rand_selection) = true;
                            end
                            foldcounter = foldcounter + 1;
                            partitions.train_indices{foldcounter} = train_idx(balanced_indices_msk);
                            partitions.test_indices{foldcounter} = fold;
                        end
                    end
                end
                    bal_partitionsall{train_modality,test_modality}{rr}=partitions; %cosmo_balance_partitions(partitions, ds_sel2);
                    ds_sel3{train_modality,test_modality}{rr} = ds_sel2;
                end

            clearvars ds_sel ds_sel2
            disp('Starting parallel process...');
            ds_sel_temp = ds_sel3{train_modality,test_modality};
            bal_part_temp = bal_partitionsall{train_modality,test_modality};

            allsamples = nan(numrepeats,241);
            allhyperplane = nan(120,241,numrepeats);
            for rr = 1:numrepeats
                ds_sel4 = ds_sel_temp{rr};
                bal_partitions = bal_part_temp{rr};
                opt=struct();
                opt.partitions=bal_partitions;
                opt.classifier=@cosmo_classify_lda3;
                opt.max_feature_count = 30000;
                opt.output = 'balanced_accuracy';
                %opt.normalization = 'zscore';

                if ~isempty(average_train_args_cell)
                    average_train_args=average_train_args_cell{1};
                    opt=cosmo_structjoin(opt, average_train_args);
                end

                measure=@cosmo_crossvalidation_measure;
    
                % Run the measure:
                % if there are not enough training samples, discard the subject:
		        % if there are not enough training samples, discard the subject:
		        if ~isempty(ds_sel4.samples)
		            num_target1_train = length(find(ds_sel4.sa.targets(bal_partitions.train_indices{1}) == 1));
		            num_target2_train = length(find(ds_sel4.sa.targets(bal_partitions.train_indices{1}) == 2));
	                if num_target1_train < ((2*trial_bin_num) + 1) || num_target2_train < ((2*trial_bin_num) + 1)
	                    disp('not enough trials! discarding subject...');
	                else
	                    nh=cosmo_interval_neighborhood(ds_sel4,'time','radius',time_radius);
				        sl_map=cosmo_searchlight_ccgp(ds_sel4,nh,measure,opt,'nproc',numCores,'progress',[]);
                        % sl_map=cosmo_searchlight_ccgp(ds_sel4,nh,measure,opt);
				        allsamples(rr,:) = sl_map.samples;
				        for t = 1:241
				            hyperplane_temp(:,t) = nanmean(cell2mat(sl_map.hyperplane{1,t}),2);
				        end
				        allhyperplane(:,:,rr) = hyperplane_temp;
	                end
		        else
		            disp('not enough trials! discarding subject...');
		        end
            end
    
	        sl_map.samples = nanmean(allsamples);
            sl_map.hyperplane = nanmean(allhyperplane,3);
	        savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Final_Resubmission_results_April16/Fig4/CCGPCrossPos/',condition,'_BasicABCD', analysistype,'_pos',num2str(posvalue1)]);
	        if ~exist(savedir,'dir'); mkdir(savedir); end
	        fprintf('Saving data for position %s pair %s combo %s %s\n', num2str(posvalue1), num2str(pair), num2str(train_modality), num2str(test_modality));
	        save(fullfile(savedir, ['/sub', num2str(s), '_pair', num2str(pair), '_combo_', ...
	            num2str(train_modality), num2str(test_modality), '_',num2str(numrepeats),'repeats_sl_map']),'sl_map','-v7.3');
	        clearvars ds_searchlight_result
	    end
    end
    end
end
end
end
% end


