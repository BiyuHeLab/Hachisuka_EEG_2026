%% Bootstrapping peak & onset times for CCGP

% Last updated June 20 2025, Ayaka Hachisuka (ahachisu@gmail.com)

clear;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
%% Load each subject:
numSubjects = 33;
subjectslist = setdiff(1:numSubjects,[22 32 33]);
bootstrapNUM = 20000;
conditionLIST = {'Normal','Masked','LSF'};
numCores = feature('numcores')-2;

trialtype = 'ALL';
time_radius = 0;        
trial_bin_num = 4;
numrepeats = 300;
n_modalities = 2;
OCCIP = 0;
run_pca = 0;
timewindow = linspace(-0.4,0.8,241);

for cc = 1:length(conditionLIST)
condition = conditionLIST{cc};

savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/resultswith200Hzdata/',condition,'_BasicABCD', trialtype]);

counter1 = 0; counter2 = 0;
for s = subjectslist
    for pair = 1:2
    if exist(fullfile([savedir,...
            '/sub',num2str(s),'_pair',num2str(pair),'_combo_11_',num2str(numrepeats),'repeats_ds_searchlight_result.mat']),'file') &&...
        exist(fullfile([savedir,...
            '/sub',num2str(s),'_pair',num2str(pair),'_combo_12_',num2str(numrepeats),'repeats_ds_searchlight_result.mat']),'file') &&...
        exist(fullfile([savedir,...
            '/sub',num2str(s),'_pair',num2str(pair),'_combo_21_',num2str(numrepeats),'repeats_ds_searchlight_result.mat']),'file') &&...
        exist(fullfile([savedir,...
            '/sub',num2str(s),'_pair',num2str(pair),'_combo_22_',num2str(numrepeats),'repeats_ds_searchlight_result.mat']),'file')
       
        for train_modality=1:n_modalities
            for test_modality=1:n_modalities
                if train_modality == test_modality
                    counter1 = counter1 + 1;
                    sl_result_group_within{s,counter1} = load(fullfile([savedir, '/sub', num2str(s), '_pair', num2str(pair), '_combo_', ...
                    num2str(train_modality), num2str(test_modality), '_',num2str(numrepeats),'repeats_ds_searchlight_result']),'ds_searchlight_result');
                else
                    counter2 = counter2 + 1;
                    sl_result_group_across{s,counter2} = load(fullfile([savedir, '/sub', num2str(s), '_pair', num2str(pair), '_combo_', ...
                    num2str(train_modality), num2str(test_modality), '_',num2str(numrepeats),'repeats_ds_searchlight_result']),'ds_searchlight_result');
                end
            end
        end
    end
    end
end

% Average within-decoding subjects:
counter = 0;

for s = 1:size(sl_result_group_within,1)
    %remove empties
    subjgrp = sl_result_group_within(s,:);
    subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
    if ~isempty(subjgrp)
        counter = counter + 1;
        submatrix = cell2mat(cellfun(@(x) x.ds_searchlight_result.samples, subjgrp, 'UniformOutput', false)');
        avg_within{1,counter}.samples = mean(submatrix,1);
        avg_within{1,counter}.a = subjgrp{1,1}.ds_searchlight_result.a;
        avg_within{1,counter}.fa = subjgrp{1,1}.ds_searchlight_result.fa;
        avg_within{1,counter}.sa = subjgrp{1,1}.ds_searchlight_result.sa;
    end
end
within_n_subj = counter;
n_subj = within_n_subj;
ds_cell = cell(n_subj,1);
for p=1:n_subj
        actual_dataset{p,1}  = avg_within{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        %round values to the 3rd decimal because "0" is a very small number
        %that jitters around 0.
        actual_dataset{p,1}.a.fdim.values{1,1} = round(avg_within{p}.a.fdim.values{1,1},2);
end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group = cosmo_stack(ds_group);

for s = 1:n_subj
    decodingacc_bysub(s,:) = avg_within{s}.samples;
end

%bootstrapping for statistics:
allow_clustering_over_time = true;
bootstrap_sig_onset = zeros(1,bootstrapNUM);
clearvars bootstrapgroup
parfor (b = 1:bootstrapNUM, numCores)
% for b = 1:bootstrapNUM
    % disp(b);
    bootstrap_indices = randi(n_subj, [1, n_subj]); %sample from 30 subjects, with replacement
    bootstrapacc = decodingacc_bysub(bootstrap_indices,:);

    bootstrapgroup = [];
    bootstrapgroup.samples = bootstrapacc;
    bootstrapgroup.sa = ds_group.sa;
    bootstrapgroup.fa = ds_group.fa;
    bootstrapgroup.a = ds_group.a;
    bootstrapgroup.a.fdim.values{1,1} = ds_group.a.fdim.values{1,1};
    bootstrapgroup.fa.center_ids = ds_group.fa.center_ids;
    bootstrapgroup.fa.time = ds_group.fa.time;

    nbrhood = cosmo_cluster_neighborhood(bootstrapgroup, 'time',allow_clustering_over_time);
    % ds_z{b} = cosmo_montecarlo_cluster_stat(bootstrapgroup,nbrhood,'h0_mean',0.5,'niter',1000,'progress',[]);
    ds_z{b} = cosmo_montecarlo_cluster_stat(bootstrapgroup,nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05,'progress',[]);
    sig = find(abs(ds_z{b}.samples)>1.96);

    if ~isempty(sig)
        bootstrap_sig_onset(b) = timewindow(sig(1)); %in ms
    else
        bootstrap_sig_onset(b) = NaN;
    end
    bootstrap_data{b} = bootstrapgroup;
end

bootstrapbycond = bootstrap_sig_onset;
bootstrap_data_all = bootstrap_data;

savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/Within/BootstrappingTwoTail',num2str(bootstrapNUM),'_LDA_',trialtype,'Trials_stats/',condition,'/'];
if ~exist(savedir,'file'); mkdir(savedir); end
save(fullfile([savedir,'bootstrap_sigval_',trialtype]),'bootstrapbycond');
save(fullfile([savedir,'bootstrap_balacc_data_',trialtype]),'bootstrap_data_all');
save(fullfile([savedir,'ds_z_',trialtype]),'ds_z');

%% across
% Average across-decoding subjects:
counter = 0;

for s = 1:size(sl_result_group_across,1)
    %remove empties
    subjgrp = sl_result_group_across(s,:);
    subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
    if ~isempty(subjgrp)
        counter = counter + 1;
        submatrix = cell2mat(cellfun(@(x) x.ds_searchlight_result.samples, subjgrp, 'UniformOutput', false)');
        avg_across{1,counter}.samples = mean(submatrix,1);
        avg_across{1,counter}.a = subjgrp{1,1}.ds_searchlight_result.a;
        avg_across{1,counter}.fa = subjgrp{1,1}.ds_searchlight_result.fa;
        avg_across{1,counter}.sa = subjgrp{1,1}.ds_searchlight_result.sa;
    end
end
across_n_subj = counter;
n_subj = across_n_subj;
ds_cell = cell(n_subj,1);
for p=1:n_subj
        actual_dataset{p,1}  = avg_across{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        %round values to the 3rd decimal because "0" is a very small number
        %that jitters around 0.
        actual_dataset{p,1}.a.fdim.values{1,1} = round(avg_across{p}.a.fdim.values{1,1},2);
end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group = cosmo_stack(ds_group);

for s = 1:n_subj
    decodingacc_bysub(s,:) = avg_across{s}.samples;
end

%bootstrapping for statistics:
allow_clustering_over_time = true;
bootstrap_sig_onset = zeros(1,bootstrapNUM);
clearvars bootstrapgroup
parfor (b = 1:bootstrapNUM, numCores)
% for b = 1:bootstrapNUM
    % disp(b);
    bootstrap_indices = randi(n_subj, [1, n_subj]); %sample from 29 subjects, with replacement
    bootstrapacc = decodingacc_bysub(bootstrap_indices,:);

    bootstrapgroup = [];
    bootstrapgroup.samples = bootstrapacc;
    bootstrapgroup.sa = ds_group.sa;
    bootstrapgroup.fa = ds_group.fa;
    bootstrapgroup.a = ds_group.a;
    bootstrapgroup.a.fdim.values{1,1} = ds_group.a.fdim.values{1,1};
    bootstrapgroup.fa.center_ids = ds_group.fa.center_ids;
    bootstrapgroup.fa.time = ds_group.fa.time;
    
    nbrhood = cosmo_cluster_neighborhood(bootstrapgroup, 'time',allow_clustering_over_time);
    % ds_z{b} = cosmo_montecarlo_cluster_stat(bootstrapgroup,nbrhood,'h0_mean',0.5,'niter',1000,'progress',[]);
    ds_z{b} = cosmo_montecarlo_cluster_stat(bootstrapgroup,nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05,'progress',[]);
    sig = find(abs(ds_z{b}.samples)>1.96);

    if ~isempty(sig)
        bootstrap_sig_onset(b) = timewindow(sig(1)); %in ms
    else
        bootstrap_sig_onset(b) = NaN;
    end
    bootstrap_data{b} = bootstrapgroup;
end

bootstrapbycond = bootstrap_sig_onset;
bootstrap_data_all = bootstrap_data;

savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/Across/BootstrappingTwoTail',num2str(bootstrapNUM),'_LDA_',trialtype,'Trials_stats/',condition,'/'];
if ~exist(savedir,'file'); mkdir(savedir); end
save(fullfile([savedir,'bootstrap_sigval_',trialtype]),'bootstrapbycond');
save(fullfile([savedir,'bootstrap_balacc_data_',trialtype]),'bootstrap_data_all');
save(fullfile([savedir,'ds_z_',trialtype]),'ds_z');
end