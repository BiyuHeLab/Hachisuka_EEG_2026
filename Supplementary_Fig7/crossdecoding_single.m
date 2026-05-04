%% MEEG CROSS-CONDITION DECODING for Normal vs. Masked, Normal vs. LSF, Masked vs. LSF

% Last updated Oct 11 2024, Ayaka Hachisuka (ahachisu@gmail.com)

% Using LDA again, Jan 16 2025

function C2F_crossdecoding_leaveoutexemp_batch(condition1,condition2,trialtype)

% cd /isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/
% cd libsvm; % change this to the directory where you put LIBSVM
% cd matlab  % go to matlab sub-directory
% make       % compile libsvm mex functions; requires a working compiler
% rmpath(pwd)   % } ensure directory is on top
% addpath(pwd)  % } of the search path
% 
% % verify it worked.
% addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
% addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');
% cosmo_check_external('libsvm'); % should not give an error

superorbasicLIST = {'Super','Basic1','Basic2'};
analysistypeLIST = {'Super','Basic','Basic'};
num_subjects = 31;
numCores = 10; %feature('numcores')-2;
numrepeats = 300;
batchsize = 10;
trial_bin_num = 1;

time_radius = 0;
run_pca = 0;
%condition1 = 'Masked'; condition2 = 'LSF';
twocondition{1} = condition1;
twocondition{2} = condition2;
%trialtype = 'Undetect';

for aa = 1:3

superorbasic = superorbasicLIST{aa}; %Basic1 or Basic2, or Super
analysistype = analysistypeLIST{aa}; %Super or Basic

if contains(superorbasic,'Super')
    categorylevel = 'Super';
elseif contains(superorbasic,'Basic')
    categorylevel = 'Basic';
elseif contains(superorbasic,'RecUnrec')
    categorylevel = 'RecUnrec';
end

if contains(superorbasic,'Super')
    index2label={'Animate','Inanimate'}; % 1=pre, 2=peri/post
elseif contains(superorbasic,'Basic1')
    index2label={'Cat','Dog'};
elseif contains (superorbasic, 'Basic2')
    index2label={'Car','Truck'};
%     index2label={'Dog','Cat','Car','Truck'}; % 1=pre, 2=peri/post
elseif contains (superorbasic, 'RecUnrec')
    index2label={'Rec1','Unrec'};
end

%% get timelock data in CoSMoMVPA format
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');

% set configuration
config=cosmo_config();
config.data_path = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_input/']);
data_path=fullfile(config.data_path);

% reset citation list
cosmo_check_external('-tic');

for s = 1:31 %24:num_subjects
    
subjNum = s;
if length(num2str(subjNum)) == 1 || ...
        (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
    subjNumStr = ['0' num2str(subjNum)];
else
    subjNumStr = num2str(subjNum);
end

% % load data
data_fn=fullfile([data_path,'sub', subjNumStr, '_dataEEG_',trialtype,'.mat']);

for c = 1:2
    condition = twocondition{c};
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
    
    if strcmp(superorbasic, 'Super')
    indices = find(data_tl.trialinfo ~= cond1 & data_tl.trialinfo ~= cond2 & data_tl.trialinfo ~= cond3 & data_tl.trialinfo ~= cond4 ...
        & data_tl.trialinfo ~= cond5 & data_tl.trialinfo ~= cond6 & data_tl.trialinfo ~= cond7 & data_tl.trialinfo ~= cond8 ...
        & data_tl.trialinfo ~= cond9 & data_tl.trialinfo ~= cond10 & data_tl.trialinfo ~= cond11 & data_tl.trialinfo ~= cond12);
    elseif strcmp(superorbasic,'Basic1')
        indices = find(data_tl.trialinfo ~= cond1 & data_tl.trialinfo ~= cond2 & data_tl.trialinfo ~= cond3 & data_tl.trialinfo ~= cond4 ...
        & data_tl.trialinfo ~= cond5 & data_tl.trialinfo ~= cond6);
    elseif strcmp(superorbasic,'Basic2')
        indices = find(data_tl.trialinfo ~= cond7 & data_tl.trialinfo ~= cond8 ...
        & data_tl.trialinfo ~= cond9 & data_tl.trialinfo ~= cond10 & data_tl.trialinfo ~= cond11 & data_tl.trialinfo ~= cond12);
    end
    
    data_tl.trialinfo(indices) = [];
    data_tl.trial(indices,:,:) = [];
    data_tl.sampleinfo(indices,:) = [];
    
    % FIND all the unique exempar numbers and re-number them, 1-12.
    [unique_values, ~, new_indices] = unique(data_tl.trialinfo);
    temptrialinfo = new_indices;
    
    data_tl.trialinfo = temptrialinfo;
    data_tlboth{c} = data_tl;
end

data_tl1 = data_tlboth{1};
data_tl2 = data_tlboth{2};
clearvars data_tlboth data_tl

% FIND all the unique exempar numbers and re-number them, 1-12 (or, 1-6 for basic)
[~, ~, new_indices] = unique(data_tl1.trialinfo);
temptrialinfo = new_indices;
data_tl1.trialinfo = temptrialinfo;

[~, ~, new_indices] = unique(data_tl2.trialinfo);
temptrialinfo = new_indices;
data_tl2.trialinfo = temptrialinfo;

% convert to cosmomvpa struct
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
ds_tl1=cosmo_meeg_dataset(data_tl1); 
ds_tl2=cosmo_meeg_dataset(data_tl2); 
ds_tl1.sa.exemptrialinfo = data_tl1.trialinfo; clearvars data_tl1
ds_tl2.sa.exemptrialinfo = data_tl2.trialinfo; clearvars data_tl2

ds_cell{1} = ds_tl1;
clearvars ds_tl1
ds_cell{2} = ds_tl2;
clearvars ds_tl2

for i = 1:2
    ds_tl = ds_cell{i};
    t = ds_tl.sa.trialinfo;
    % TARGETS
    if strcmp(superorbasic,'Super')
        ds_tl.sa.targets( t == 1 | t == 2 | t == 3 | t == 4 | t == 5 | t == 6) = 1; % animate
        ds_tl.sa.targets(t == 7 | t == 8 | t == 9 | t == 10 | t == 11 | t == 12) = 2; % inanimate
        ds_tl.sa.targets = ds_tl.sa.targets';
    elseif contains(superorbasic,'Basic')
        ds_tl.sa.targets(t == 1 | t == 2 | t == 3) = 1; % Dog or Car
        ds_tl.sa.targets(t == 4 | t == 5 | t == 6) = 2; % Cat or Truck
        ds_tl.sa.targets = ds_tl.sa.targets';
    end

    % MODALITY
    ds_tl.sa.modality = ones(size(ds_tl.samples,1),1) + (i-1); %1 for cond1, 2 for cond2.

    ds_cell{i} = ds_tl;
end

clearvars ds_tl

% Combine datasets from two conditions:
[idxs,ds_intersect_cell]=cosmo_mask_dim_intersect(ds_cell);
clearvars ds_cell
ds_all=cosmo_stack(ds_intersect_cell,1);
clearvars ds_intersect_cell
ds_all.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_all.sa.targets));

% just to check everything is ok
cosmo_check_dataset(ds_all);

%% do all combinations for training and test modalities

n_modalities=2;
totalrepeat = numrepeats;

partitionsall = cell(n_modalities,n_modalities);
ds_sel3 = cell(n_modalities, n_modalities);
sl_map_samples = cell(n_modalities,n_modalities);

for train_modality=1:n_modalities
    for test_modality=1:n_modalities
        sl_map_samples{train_modality,test_modality} = nan(numrepeats,241);
        for batch = 1:(numrepeats/batchsize)
            % Start and end indices for this batch
            startIdx = (batch - 1) * batchsize + 1;
            endIdx = batch * batchsize;
            
            if train_modality~=test_modality
                
                msk = cosmo_match(ds_all.sa.modality,[train_modality test_modality]);
                ds_sel=cosmo_slice(ds_all,msk);
                
                partitionsall{train_modality,test_modality} = cell((numrepeats/batchsize),1);
                ds_sel3{train_modality,test_modality} = cell((numrepeats/batchsize),1);

                for repeat = startIdx:endIdx
                %balance exemplar count:
                [unique_values, ~, indices] = unique(ds_sel.sa.exemptrialinfo);
                exemp_counts = accumarray(indices(:), 1);  % how many times does each exemplar come up
                min_exemp_count = min(exemp_counts);
                selected_trials = false(size(ds_sel.sa.exemptrialinfo)); %logical index mask

                % Loop over each unique exemplar and reduce trial count to
                % minimum exemplar count:
                for i = unique_values'
                    exemp_ind = find(ds_sel.sa.exemptrialinfo == unique_values(i));
                    rand_selection = randsample(exemp_ind, min_exemp_count); %sample without replacement to minimize value 
                    selected_trials(rand_selection) = true;
                end
                
                balanced_data = ds_sel.sa.exemptrialinfo(selected_trials);

                %remove all the other trials
                balanced_indices_msk = selected_trials;

                ds_sel2 = cosmo_slice(ds_sel,balanced_indices_msk);
                ds_sel2.sa.chunks = ds_sel2.sa.modality; %(1:size(ds_sel.samples,1))';
                notenoughinmodality = 0;
                if ~all(ismember([1, 2], ds_sel2.sa.modality))
                    disp('Not enough trials in each modality...');
                    notenoughinmodality = 1;
                    break;
                end

                allindex = 1:size(ds_sel2.samples,1); %vector of all trial indices
 
                % Partition the balanced data:
                partitions = cosmo_nchoosek_partitioner(ds_sel2,1,'modality',test_modality);
                partitionsall{train_modality,test_modality}{repeat} = partitions;
                ds_sel3{train_modality,test_modality}{repeat} = ds_sel2;

                % if there are not enough training samples, discard the subject:
                num_target1_train = length(find(ds_sel2.sa.targets(allindex(partitions.train_indices{1})) == 1));
                num_target2_train = length(find(ds_sel2.sa.targets(allindex(partitions.train_indices{1})) == 2));
                %you need at least 2 samples per class in training:
                notenoughtrial{train_modality,train_modality} = 0;
                if num_target1_train < ((2*trial_bin_num) + 1) || num_target2_train < ((2*trial_bin_num) + 1)
                %     repeat = repeat-1;
                %     totalrepeat = totalrepeat+1;
                % end
                % 
                % if totalrepeat > 20
                    disp('Not enough trials...');
                    notenoughtrial{train_modality,train_modality} = 1;
                    break;
                end
            end

clearvars ds_sel ds_sel2

% Run PARFOR loop across batches:
        if notenoughtrial{train_modality,train_modality} == 0 && notenoughinmodality == 0
            sl_map = cell(numrepeats,1);
            partitions = cell(numrepeats,1);
            tempsamples = zeros(numrepeats,241);
            parfor (rr = startIdx:endIdx,numCores)
            % for rr = startIdx:endIdx
                ds_sel4 = ds_sel3{train_modality,test_modality}{rr};
                partitions{rr}.test_indices{1} =  partitionsall{train_modality,test_modality}{rr}.test_indices{1};
                partitions{rr}.train_indices{1} =  partitionsall{train_modality,test_modality}{rr}.train_indices{1};

                bal_partitions{rr}=cosmo_balance_partitions(partitions{rr}, ds_sel4);
                % in addition give a label to each trial
                ds_sel4.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_sel4.sa.targets));
                
                % just to check everything is ok
                cosmo_check_dataset(ds_sel4);
                
                average_train_args_cell={{'average_train_count',trial_bin_num,...            %  { 2
                                                'average_train_resamplings',1}}; %  { 2
                n_average_train_args=numel(average_train_args_cell);
                
                chantypes = {'eeg'};
                
                nchantypes=length(chantypes);
                
                npartitions=numel(bal_partitions{rr});
                fprintf('There are %d partitions\n', numel(bal_partitions{rr}.train_indices));
                fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                                        bal_partitions{rr}.train_indices)));
                fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                                        bal_partitions{rr}.test_indices)));
                measure=@cosmo_crossvalidation_measure_orig;
                measure_args=struct();
                % measure_args.classifier=@cosmo_classify_libsvm;
                measure_args.classifier=@cosmo_classify_lda3
                measure_args.partitions=bal_partitions{rr};
                measure_args.output = 'balanced_accuracy';

                if run_pca == 1
                    measure_args.pca_explained_ratio = 0.95;
                end
                % add the options to average samples to the measure arguments.
                % (if no averaging is desired, this step can be left out.)
                if ~isempty(average_train_args_cell)
                    average_train_args=average_train_args_cell{1};
                    measure_args=cosmo_structjoin(measure_args, average_train_args);
                end
                                      
                nbrhood=cosmo_interval_neighborhood(ds_sel4,'time',...
                                                        'radius',time_radius);
                sl_map{rr}=cosmo_searchlight(ds_sel4,nbrhood,measure,measure_args);
                tempsamples(rr,:) = sl_map{rr}.samples; %save every repeat accuracy
            end
            sl_map_samples{train_modality,test_modality}(startIdx:endIdx,:) = tempsamples(startIdx:endIdx,:);
        end
        end
        end
    end
end

if exist('sl_map','var')
for train_modality=1:n_modalities
    for test_modality=1:n_modalities
        if train_modality ~= test_modality
            ds_searchlight_result = sl_map{1};
            ds_searchlight_result.samples = nanmean(sl_map_samples{train_modality,test_modality});
            savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'/']);
            if ~exist(savedir); mkdir(savedir); end
            save([savedir,twocondition{1},'_',twocondition{2},'sub',num2str(s), '_',...
                num2str(train_modality), num2str(test_modality), '_',num2str(numrepeats),'repeats_ds_searchlight_result'],'ds_searchlight_result','-v7.3');
        end
    end
end
end
end
end
end
