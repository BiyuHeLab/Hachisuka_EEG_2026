%% Plotting PAS vs. Distance to Bound (LDA)

% This script plots the following (with outputs from
% C2F_contentdecoding_distancetobound_bytrial_withrepeatv2.m)
% 1) Raw distance values, z-scored, baseline-corrected, and smoothed using
% the Savitzky-Golay filter with a 100ms window.The plot averages
% single-trial level distances by PAS rating, plotting for traces total for
% each of the ratings.
% 2) Decoding accuracy outputs from this leave-one-trial out schematic.
% For bootstrapping results of (2), see C2F_dist2bound_acc_bootstrap.m
% 3) An option to plot the decoding accuracy outputs by PAS rating (not
% included in the MS).

% Last updated July 1 2025, Ayaka Hachisuka (ahachisu@gmail.com)

clear; 
num_subjects = 31;
subjectslist = setdiff(1:num_subjects,[22 32 33]);

whichcondA = 1:3:36;
whichcondB = 2:3:36;
whichcondC = 3:3:36;
whichcondA1 = 1:3:18;
whichcondB1 = 2:3:18;
whichcondC1 = 3:3:18;
whichcondA2 = 19:3:36;
whichcondB2 = 20:3:36;
whichcondC2 = 21:3:36;

cond1 = 'AllCond'; %AllCond_50Hz';
cond = 'AllCond';
trialtype = 'ALL';

load(['/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/Behavioral/C2F_BehavData_30Subs.mat']);
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');

% set configuration
config=cosmo_config();
config.data_path = fullfile('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_inputv2/');
data_path=fullfile(config.data_path);

%% Load distance to bound

time_radius = 0;        
trial_bin_num = 1;
trial_resample = 1;
nminval = 1;

config.data_path = fullfile('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/PASrating1/');
data_path=fullfile(config.data_path);

%% Super
load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/Behavioral/trialinfo_alltrials_30subsV2.mat','trialinfo');
save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/Super/noPCA_35Hz_', num2str(time_radius), ...
    'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
counter = 0;
trialsofinterest = 1:36;
alltrialcount = 0;
for s = 1:length(subjectslist)
    
    trials = find(ismember(trialinfo{s}, trialsofinterest));
    clearvars pasresp pas0ind pas1ind pas2ind pas3ind pasresp3 pasresp2 condbytrial supercorr basiccorr
    
    counter = counter + 1;
   
    subjNum = subjectslist(s);
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end

    load([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',cond,'_mvpaoutput'],'sl_map2');

    data_fn=fullfile([data_path,'sub', num2str(subjNumStr), '_dataEEG_',trialtype,'.mat']);
    data_tl1=load(data_fn);

    sl_map = sl_map2;
    sl_map.origtrialinfo = data_tl1.data3.origtrialinfo;

    trialnumbers = NaN(size(sl_map.dist,1), 1);                
    valid_rows = any(~isnan(sl_map.dist), 2); % rows without NaN
    trialnumbers(valid_rows) = find(valid_rows); % Add index is the row is valid.

    origtrialinfo = sl_map.origtrialinfo;
    pasresp = origtrialinfo(:,19);
    condbytrial = origtrialinfo(:,9);
    supercorr = origtrialinfo(:,15);
    basiccorr = origtrialinfo(:,17);
    exempbytrial = origtrialinfo(:,12);
    basiccorr(basiccorr<0) = 0;
    exemptrials = origtrialinfo(:,12);
    %remove trials where PAS is NaN (no response):
    nantrials = find(isnan(pasresp));
    if ~isempty(nantrials)
        trialnumbers(nantrials) = [];
        pasresp(nantrials) = []; 
        condbytrial(nantrials) = [];
        sl_map.dist(nantrials,:) = [];
        sl_map.samples(nantrials,:) = [];
        exemptrials(nantrials,:) = [];
        supercorr(nantrials,:) = [];
        basiccorr(nantrials,:) = [];
        exempbytrial(nantrials,:) = [];
    end

    NormalInd = find(condbytrial == 1);
    MaskedInd = find(condbytrial == 2);
    LSFInd = find(condbytrial == 3);

    % SIGNED DISTANCE
    trialinfo2 = nan(size(exemptrials));
    trialinfo2(exemptrials <= 6) = 1; % animal trials
    trialinfo2(exemptrials > 6) = 2; % vehicle trials
    trialinfo2_expanded = repmat(trialinfo2, 1, 241);
    % if it's an inanimate trial, multiply by -1. This keeps the negative
    % sign of misclassified distances.
    sl_map.dist(trialinfo2_expanded == 2) = -1 * sl_map.dist(trialinfo2_expanded == 2);

     % correct baseline
     % baseline = nanmean(sl_map.dist(:,1:81),2);
     % sl_map.dist = sl_map.dist - baseline;

    ecounter = 0;
    for exemp = trialsofinterest
        ecounter = ecounter + 1;
        exemptrial = find(ismember(exempbytrial, exemp));
        distbyexemp(ecounter,:) = nanmean(sl_map.dist(exemptrial,:));
    end

    pas0ind = find(pasresp == 0);
    pas1ind = find(pasresp == 1);
    pas2ind = find(pasresp == 2);
    pas3ind = find(pasresp == 3);

    corrind = find(supercorr == 1);
    incorrind = find(supercorr == 0);

    keep_trials = 1:length(condbytrial);

    alldist_lmm{s} = sl_map.dist;
    allpas_lmm{s} = pasresp;
    condtype_lmm{s} = condbytrial(keep_trials, :);
    allexemp_lmm{s} = exempbytrial(keep_trials, :);
    trialnum_lmm{s} = 1:length(pasresp);
    allsupercorrect_lmm{s} = supercorr(keep_trials, :);
    allbasiccorrect_lmm{s} = basiccorr(keep_trials, :);

    dist_for_stats{1}{s}(1,:) = nanmean((sl_map.dist(pas0ind,:)),1);
    dist_for_stats{1}{s}(2,:) = nanmean((sl_map.dist(pas1ind,:)),1);
    dist_for_stats{1}{s}(3,:) = nanmean((sl_map.dist(pas2ind,:)),1);
    dist_for_stats{1}{s}(4,:) = nanmean((sl_map.dist(pas3ind,:)),1);

    acc_for_stats{1}{s}(1,:) = nanmean((sl_map.samples(pas0ind,:)),1);
    acc_for_stats{1}{s}(2,:) = nanmean((sl_map.samples(pas1ind,:)),1);
    acc_for_stats{1}{s}(3,:) = nanmean((sl_map.samples(pas2ind,:)),1);
    acc_for_stats{1}{s}(4,:) = nanmean((sl_map.samples(pas3ind,:)),1);

    alldist_basic1(:,:,counter) = distbyexemp; % (distbyexemp - min(distbyexemp,[],2)) ./ (max(distbyexemp,[],2) - min(distbyexemp,[],2)); %distbyexemp;
    alldist0(counter,:) = nanmean((sl_map.dist(pas0ind,:)),1);
    alldist1(counter,:) = nanmean((sl_map.dist(pas1ind,:)),1);
    alldist2(counter,:) = nanmean((sl_map.dist(pas2ind,:)),1);
    alldist3(counter,:) = nanmean((sl_map.dist(pas3ind,:)),1);

    alldist_corr(counter,:) = nanmean((sl_map.dist(corrind,:)),1);
    alldist_incorr(counter,:) = nanmean((sl_map.dist(incorrind,:)),1);
    allacc1(counter,:) = nanmean((sl_map.samples));

end

figure; hold on
plot(nanmean(alldist_corr),'color','b');
plot(nanmean(alldist_incorr),'color','r');
% plot(nanmean(alldist2));
% plot(nanmean(alldist3),'r');

%%
lmm_super = cell(241,1);
for t = 1:241
    pas_lmm = [];
    supercorr_lmm = [];
    basiccorr_lmm = [];
    dist_lmm = [];
    exemp_lmm = [];
    cond_lmm = [];
    trial_lmm = [];
    subj_lmm = [];
    for s = 1:length(subjectslist)
        supercorr_lmm = [supercorr_lmm; allsupercorrect_lmm{s}];
        basiccorr_lmm = [basiccorr_lmm; allbasiccorrect_lmm{s}];

        pas_lmm = [pas_lmm; allpas_lmm{s}];        
        dist_lmm = [dist_lmm; alldist_lmm{s}(:,t)];  
        exemp_lmm = [exemp_lmm; allexemp_lmm{s}];
        cond_lmm = [cond_lmm; condtype_lmm{s}];       
        trial_lmm = [trial_lmm; trialnum_lmm{s}'];       
        subj_lmm = [subj_lmm; repmat(s,length(allpas_lmm{s}),1)];       
    end
    pas_lmm(pas_lmm == 3) = 4;
    pas_lmm(pas_lmm == 2) = 3;
    pas_lmm(pas_lmm == 1) = 2;
    pas_lmm(pas_lmm == 0) = 1;

    lmm_super_table{t} = table(categorical(supercorr_lmm), categorical(basiccorr_lmm), ...
        categorical(pas_lmm), dist_lmm, categorical(cond_lmm), categorical(exemp_lmm), trial_lmm, subj_lmm, ...
        'VariableNames', {'SuperCorr','BasicCorr','PAS', 'Dist', 'Cond', 'Exemp','Trial', 'Subject'});
end

%% Basic 1
clearvars distbyexemp alldist_lmm
trialsofinterest = 1:18;
load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/Behavioral/trialinfo_alltrials_30subsV2.mat','trialinfo');
save_dir1=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/Basic1/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
counter = 0;
for s = 1:length(subjectslist)
    clearvars pasresp pas0ind pas1ind pas2ind pas3ind pasresp3 pasresp2 condbytrial supercorr basiccorr trialnumbers
    counter = counter + 1;

    subjNum = subjectslist(s);
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end
    load([save_dir1,'sub',num2str(subjNumStr),'_dataEEG_',cond,'_mvpaoutput'],'sl_map2');
    data_fn=fullfile([data_path,'sub', num2str(subjNumStr), '_dataEEG_',trialtype,'.mat']);
    data_tl1=load(data_fn);

    sl_map = sl_map2;
    sl_map.origtrialinfo = data_tl1.data3.origtrialinfo;

    if s == 30
        A = sl_map.dist;
        [max_val, linear_idx] = max(abs(A(:)));
        [row_idx, col_idx] = ind2sub(size(A), linear_idx);

        sl_map.dist(row_idx,:) = [];
        sl_map.samples(row_idx,:) = [];
        sl_map.origtrialinfo(row_idx,:) = [];
    end

    trialnumbers = NaN(size(sl_map.dist,1), 1);                
    valid_rows = any(~isnan(sl_map.dist), 2); % rows without NaN
    trialnumbers(valid_rows) = find(valid_rows); % Add index is the row is valid.

    origtrialinfo = sl_map.origtrialinfo;
    pasresp = origtrialinfo(:,19);
    condbytrial = origtrialinfo(:,9);
    supercorr = origtrialinfo(:,15);
    basiccorr = origtrialinfo(:,17);
    exempbytrial = origtrialinfo(:,12);
    basiccorr(basiccorr<0) = 0;
    exemptrials = origtrialinfo(:,12);
    %remove trials where PAS is NaN (no response):
    nantrials = find(isnan(pasresp));
    if ~isempty(nantrials)
        trialnumbers(nantrials) = [];
        pasresp(nantrials) = []; 
        condbytrial(nantrials) = [];
        sl_map.dist(nantrials,:) = [];
        sl_map.samples(nantrials,:) = [];
        exemptrials(nantrials,:) = [];
        supercorr(nantrials,:) = [];
        basiccorr(nantrials,:) = [];
        exempbytrial(nantrials,:) = [];
    end

    pas0ind = find(pasresp == 0);
    pas1ind = find(pasresp == 1);
    pas2ind = find(pasresp == 2);
    pas3ind = find(pasresp == 3);

    corrind = find(basiccorr == 1);
    incorrind = find(basiccorr == 0);

    NormalInd = find(condbytrial == 1);
    MaskedInd = find(condbytrial == 2);
    LSFInd = find(condbytrial == 3);

    % SIGNED DISTANCE
    trialinfo2 = nan(size(exemptrials));
    trialinfo2(exemptrials <= 3) = 1; % dog trials
    trialinfo2(exemptrials > 3) = 2; % cat trials
    trialinfo2_expanded = repmat(trialinfo2, 1, 241);
    % if it's an inanimate trial, multiply by -1. This keeps the negative
    % sign of misclassified distances.
    sl_map.dist(trialinfo2_expanded == 2) = -1 * sl_map.dist(trialinfo2_expanded == 2);

    % correct baseline
    % baseline = nanmean(sl_map.dist(:,1:81),2);
    % sl_map.dist = sl_map.dist - baseline;

    ecounter = 0;
    for exemp = trialsofinterest
        ecounter = ecounter + 1;
        exemptrial = find(ismember(exempbytrial, exemp));
        distbyexemp(ecounter,:) = nanmean(sl_map.dist(exemptrial,:));
    end

    % rearrage dist by exemplar:
    reorganizedtrials = [1:3:18,2:3:18,3:3:18];
    distbyexemp = distbyexemp(reorganizedtrials,:);
    distbyexemp_basic{1}(:,:,counter) = distbyexemp;
    
    pas0ind = find(pasresp == 0);
    pas1ind = find(pasresp == 1);
    pas2ind = find(pasresp == 2);
    pas3ind = find(pasresp == 3);
    keep_trials = 1:length(condbytrial);
    alldist_lmm{s} = sl_map.dist;
    allpas_lmm{s} = pasresp;
    condtype_lmm{s} = condbytrial(keep_trials, :);
    allexemp_lmm{s} = exempbytrial(keep_trials, :);
    trialnum_lmm{s} = 1:length(pasresp);
    allsupercorrect_lmm{s} = supercorr(keep_trials, :);
    allbasiccorrect_lmm{s} = basiccorr(keep_trials, :);
    
    dist_for_stats{2}{s}(1,:) = nanmean((sl_map.dist(pas0ind,:)),1);
    dist_for_stats{2}{s}(2,:) = nanmean((sl_map.dist(pas1ind,:)),1);
    dist_for_stats{2}{s}(3,:) = nanmean((sl_map.dist(pas2ind,:)),1);
    dist_for_stats{2}{s}(4,:) = nanmean((sl_map.dist(pas3ind,:)),1);

    acc_for_stats{2}{s}(1,:) = nanmean((sl_map.samples(pas0ind,:)),1);
    acc_for_stats{2}{s}(2,:) = nanmean((sl_map.samples(pas1ind,:)),1);
    acc_for_stats{2}{s}(3,:) = nanmean((sl_map.samples(pas2ind,:)),1);
    acc_for_stats{2}{s}(4,:) = nanmean((sl_map.samples(pas3ind,:)),1);

    alldist0a(counter,:) = nanmean((sl_map.dist(pas0ind,:)),1);
    alldist1a(counter,:) = nanmean((sl_map.dist(pas1ind,:)),1);
    alldist2a(counter,:) = nanmean((sl_map.dist(pas2ind,:)),1);
    alldist3a(counter,:) = nanmean((sl_map.dist(pas3ind,:)),1);
    allacc1(counter,:) = nanmean((sl_map.samples));

    alldist_corr(counter,:) = nanmean((sl_map.dist(corrind,:)),1);
    alldist_incorr(counter,:) = nanmean((sl_map.dist(incorrind,:)),1);

end

figure; hold on
plot(nanmean([alldist_corr]));
plot(nanmean([alldist_incorr]));
% plot(nanmean([alldist2a]));
% plot(nanmean([alldist3a]),'r');

%%
lmm_basic1 = cell(241,1);
for t = 1:241
    pas_lmm = [];
    supercorr_lmm = [];
    basiccorr_lmm = [];
    dist_lmm = [];
    exemp_lmm = [];
    cond_lmm = [];
    trial_lmm = [];
    subj_lmm = [];
    for s = 1:length(subjectslist)
        supercorr_lmm = [supercorr_lmm; allsupercorrect_lmm{s}];
        basiccorr_lmm = [basiccorr_lmm; allbasiccorrect_lmm{s}];

        pas_lmm = [pas_lmm; allpas_lmm{s}];        
        dist_lmm = [dist_lmm; alldist_lmm{s}(:,t)];  
        exemp_lmm = [exemp_lmm; allexemp_lmm{s}];
        cond_lmm = [cond_lmm; condtype_lmm{s}];       
        trial_lmm = [trial_lmm; trialnum_lmm{s}'];       
        subj_lmm = [subj_lmm; repmat(s,length(allpas_lmm{s}),1)];       
    end
    pas_lmm(pas_lmm == 3) = 4;
    pas_lmm(pas_lmm == 2) = 3;
    pas_lmm(pas_lmm == 1) = 2;
    pas_lmm(pas_lmm == 0) = 1;

    lmm_basic1_table{t} = table(categorical(supercorr_lmm), categorical(basiccorr_lmm), ...
        categorical(pas_lmm), dist_lmm, categorical(cond_lmm), categorical(exemp_lmm), trial_lmm, subj_lmm, ...
        'VariableNames', {'SuperCorr','BasicCorr','PAS', 'Dist', 'Cond', 'Exemp','Trial', 'Subject'});
end

%% Basic 2
clearvars distbyexemp alldist_lmm
trialsofinterest = 19:36;
load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/Behavioral/trialinfo_alltrials_30subsV2.mat','trialinfo');
save_dir2=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/Basic2/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
counter = 0;
for s = 1:length(subjectslist)
    clearvars pasresp pas0ind pas1ind pas2ind pas3ind pasresp3 pasresp2 condbytrial supercorr basiccorr
    counter = counter + 1;

    subjNum = subjectslist(s);
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end
    load([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',cond,'_mvpaoutput'],'sl_map2');
    data_fn=fullfile([data_path,'sub', num2str(subjNumStr), '_dataEEG_',trialtype,'.mat']);
    data_tl1=load(data_fn);

    sl_map = sl_map2;
    sl_map.origtrialinfo = data_tl1.data3.origtrialinfo;
    
    trialnumbers = NaN(size(sl_map.dist,1), 1);                
    valid_rows = any(~isnan(sl_map.dist), 2); % rows without NaN
    trialnumbers(valid_rows) = find(valid_rows); % Add index is the row is valid.

    origtrialinfo = sl_map.origtrialinfo;
    pasresp = origtrialinfo(:,19);
    condbytrial = origtrialinfo(:,9);
    supercorr = origtrialinfo(:,15);
    basiccorr = origtrialinfo(:,17);
    exempbytrial = origtrialinfo(:,12);
    basiccorr(basiccorr<0) = 0;
    exemptrials = origtrialinfo(:,12);
    %remove trials where PAS is NaN (no response):
    nantrials = find(isnan(pasresp));
    if ~isempty(nantrials)
        trialnumbers(nantrials) = [];
        pasresp(nantrials) = []; 
        condbytrial(nantrials) = [];
        sl_map.dist(nantrials,:) = [];
        sl_map.samples(nantrials,:) = [];
        exemptrials(nantrials,:) = [];
        supercorr(nantrials,:) = [];
        basiccorr(nantrials,:) = [];
        exempbytrial(nantrials,:) = [];
    end

    pas0ind = find(pasresp == 0);
    pas1ind = find(pasresp == 1);
    pas2ind = find(pasresp == 2);
    pas3ind = find(pasresp == 3);

    NormalInd = find(condbytrial == 1);
    MaskedInd = find(condbytrial == 2);
    LSFInd = find(condbytrial == 3);

    % SIGNED DISTANCE
    trialinfo2 = nan(size(exemptrials));
    trialinfo2(exemptrials <= 9) = 1; % dog trials
    trialinfo2(exemptrials > 9) = 2; % cat trials
    trialinfo2_expanded = repmat(trialinfo2, 1, 241);
    % if it's an inanimate trial, multiply by -1. This keeps the negative
    % sign of misclassified distances.
    sl_map.dist(trialinfo2_expanded == 2) = -1 * sl_map.dist(trialinfo2_expanded == 2);

     % correct baseline
    % baseline = nanmean(sl_map.dist(:,1:81),2);
    % sl_map.dist = sl_map.dist - baseline;

    ecounter = 0;
    for exemp = trialsofinterest
        ecounter = ecounter + 1;
        exemptrial = find(ismember(exempbytrial, exemp));
        distbyexemp(ecounter,:) = nanmean(sl_map.dist(exemptrial,:));
    end
    % meanpasall2 = meanpasall(:,[7:12,19:24,31:36]);

    % rearrage dist by exemplar:
    reorganizedtrials = [1:3:18,2:3:18,3:3:18];
    distbyexemp = distbyexemp(reorganizedtrials,:);
    distbyexemp_basic{2}(:,:,counter) = distbyexemp;

    pas0ind = find(pasresp == 0);
    pas1ind = find(pasresp == 1);
    pas2ind = find(pasresp == 2);
    pas3ind = find(pasresp == 3);

    corrind = find(basiccorr == 1);
    incorrind = find(basiccorr == 0);

    keep_trials = 1:length(condbytrial);
    alldist_lmm{s} = sl_map.dist;
    allpas_lmm{s} = pasresp;
    condtype_lmm{s} = condbytrial(keep_trials, :);
    trialnum_lmm{s} = 1:length(pasresp);
    allexemp_lmm{s} = exempbytrial(keep_trials, :);
    allsupercorrect_lmm{s} = supercorr(keep_trials, :);
    allbasiccorrect_lmm{s} = basiccorr(keep_trials, :);
    
    dist_for_stats{3}{s}(1,:) = nanmean((sl_map.dist(pas0ind,:)),1);
    dist_for_stats{3}{s}(2,:) = nanmean((sl_map.dist(pas1ind,:)),1);
    dist_for_stats{3}{s}(3,:) = nanmean((sl_map.dist(pas2ind,:)),1);
    dist_for_stats{3}{s}(4,:) = nanmean((sl_map.dist(pas3ind,:)),1);

    acc_for_stats{3}{s}(1,:) = nanmean((sl_map.samples(pas0ind,:)),1);
    acc_for_stats{3}{s}(2,:) = nanmean((sl_map.samples(pas1ind,:)),1);
    acc_for_stats{3}{s}(3,:) = nanmean((sl_map.samples(pas2ind,:)),1);
    acc_for_stats{3}{s}(4,:) = nanmean((sl_map.samples(pas3ind,:)),1);

    alldist_basic2(:,:,counter) = distbyexemp; % (distbyexemp - min(distbyexemp,[],2)) ./ (max(distbyexemp,[],2) - min(distbyexemp,[],2)); %distbyexemp;
    alldist0b(counter,:) = nanmean((sl_map.dist(pas0ind,:)),1);
    alldist1b(counter,:) = nanmean((sl_map.dist(pas1ind,:)),1);
    alldist2b(counter,:) = nanmean((sl_map.dist(pas2ind,:)),1);
    alldist3b(counter,:) = nanmean((sl_map.dist(pas3ind,:)),1);
    allacc2(counter,:) = nanmean((sl_map.samples));

    alldist_corr(counter,:) = nanmean((sl_map.dist(corrind,:)),1);
    alldist_incorr(counter,:) = nanmean((sl_map.dist(incorrind,:)),1);
end

figure; hold on
plot(movmean(nanmean(alldist_corr),1));
plot(movmean(nanmean(alldist_incorr),1));
% plot(movmean(nanmean(alldist2b),1));
% plot(movmean(nanmean(alldist3b),1),'r');
%%
lmm_basic2 = cell(241,1);
for t = 1:241
    pas_lmm = [];
    supercorr_lmm = [];
    basiccorr_lmm = [];
    dist_lmm = [];
    exemp_lmm = [];
    cond_lmm = [];
    trial_lmm = [];
    subj_lmm = [];
    for s = 1:length(subjectslist)
        supercorr_lmm = [supercorr_lmm; allsupercorrect_lmm{s}];
        basiccorr_lmm = [basiccorr_lmm; allbasiccorrect_lmm{s}];

        pas_lmm = [pas_lmm; allpas_lmm{s}];        
        dist_lmm = [dist_lmm; alldist_lmm{s}(:,t)];  
        exemp_lmm = [exemp_lmm; allexemp_lmm{s}];
        cond_lmm = [cond_lmm; condtype_lmm{s}];       
        trial_lmm = [trial_lmm; trialnum_lmm{s}'];       
        subj_lmm = [subj_lmm; repmat(s,length(allpas_lmm{s}),1)];     
    end
    pas_lmm(pas_lmm == 3) = 4;
    pas_lmm(pas_lmm == 2) = 3;
    pas_lmm(pas_lmm == 1) = 2;
    pas_lmm(pas_lmm == 0) = 1;

    lmm_basic2_table{t} = table(categorical(supercorr_lmm), categorical(basiccorr_lmm), ...
        categorical(pas_lmm), dist_lmm, categorical(cond_lmm), categorical(exemp_lmm), trial_lmm, subj_lmm, ...
        'VariableNames', {'SuperCorr','BasicCorr','PAS', 'Dist', 'Cond', 'Exemp','Trial', 'Subject'});
end
%% Re-organizing variables and saving them for CLMM:
for t = 1:241
lmm_super_table{t}.PAS = double(lmm_super_table{t}.PAS);
lmm_basic1_table{t}.PAS = double(lmm_basic1_table{t}.PAS);
lmm_basic2_table{t}.PAS = double(lmm_basic2_table{t}.PAS);
end

lmm_super_crosstime = table();
for t = 1:241
    temp = lmm_super_table{t};
    temp.timebin = repmat(t, height(temp), 1);
    lmm_super_crosstime = [lmm_super_crosstime; temp];
end

lmm_basic_crosstime = table();
for t = 1:241
    temp = [lmm_basic1_table{t}; lmm_basic2_table{t}];
    temp.timebin = repmat(t, height(temp), 1);
    lmm_basic_crosstime = [lmm_basic_crosstime; temp];
end

lmm_basic_crosstime.NewTrial = zeros(height(lmm_basic_crosstime), 1);
subjects = unique(lmm_basic_crosstime.Subject);

for i = 1:length(subjects)
    subj_idx = lmm_basic_crosstime.Subject == subjects(i);
    timebins = unique(lmm_basic_crosstime.timebin(subj_idx));
    
    for j = 1:length(timebins)
        idx = find(subj_idx & lmm_basic_crosstime.timebin == timebins(j));
        trials = lmm_basic_crosstime.Trial(idx);
        
        % Detect start of blocks (when trial number resets or decreases)
        resets = [true; diff(trials) < 0];  % start of each block
        block_ids = cumsum(resets);        % block counter
        
        max_trial = max(trials);
        lmm_basic_crosstime.NewTrial(idx) = trials + (block_ids - 1) * max_trial;
    end
end

lmm_basic_crosstime.Trial = [];
lmm_basic_crosstime.Trial = lmm_basic_crosstime.NewTrial;

% save data

writetable(lmm_super_crosstime, ['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/lmm_super.csv']);
writetable(lmm_basic_crosstime, ['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/lmm_basic.csv']);

%% Average Basic 1 & 2:
basic1dist = dist_for_stats{2};
basic2dist = dist_for_stats{3};

basic1acc = acc_for_stats{2};
basic2acc = acc_for_stats{3};

avg_dist_data = cell(length(subjectslist),1);
avg_acc_data = cell(length(subjectslist),1);
for i = 1:length(subjectslist)
    avg_dist_data{i} = (basic1dist{i} + basic2dist{i}) / 2;
    avg_acc_data{i} = (basic1acc{i} + basic2acc{i}) / 2;
end

dist_for_stats{2} = [];
acc_for_stats{2} = [];

dist_for_stats{2} = avg_dist_data;
acc_for_stats{2} = avg_acc_data;

%% Decoding analysis: all trials


plotcolor{1} = [50 130 246]./255;
plotcolor{2} = [240 134 80]./255;

y_pos(1) = 0.475;
y_pos(2) = 0.485;

timewindow = linspace(-0.4,0.8,241); %-0.4:0.01:0.8;
figure; set(gcf,'Color','w');
hold on;

for i = 1:2


clearvars actual_dataset ds_group
for p=1:length(subjectslist)
    if i == 1
        actual_dataset{p,1}.samples = nanmean(acc_for_stats{i}{p});
    else
        actual_dataset{p,1}.samples = (nanmean(acc_for_stats{2}{p}) + nanmean(acc_for_stats{3}{p}))./2;
    end
    actual_dataset{p,1}.sa.targets = 1;
    actual_dataset{p,1}.sa.chunks = p;
    actual_dataset{p,1}.fa.time = (1:241);
    actual_dataset{p,1}.fa.center_ids = (1:241);
    actual_dataset{p,1}.a.meeg.samples_field = 'trial';
    actual_dataset{p,1}.a.fdim.labels = {'time'};
    actual_dataset{p,1}.a.fdim.values{1,1} = linspace(-0.4,0.8,241);

    data_for_bf_clmm{i}(p,:) = actual_dataset{p}.samples;
end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group = cosmo_stack(ds_group);

allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group, 'time',allow_clustering_over_time);
ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);
sig_within = find(abs(ds_z.samples)>1.96);
sigmask = abs(ds_z.samples)>1.96;

clearvars temp
plottime = timewindow;

for s = 1:length(subjectslist)
    temp(s,:) = actual_dataset{s,1}.samples;
end

meansub = nanmean(temp);
sem_sub = nanstd(temp)/sqrt(length(subjectslist));
x = timewindow;
curve1 = meansub + sem_sub;
curve2 = meansub - sem_sub;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

h = fill(x2, inBetween, plotcolor{i});
set(h, 'EdgeColor', 'none'); % Remove black border
alpha(0.1);

time_values = timewindow;
plot1 = plot(plottime,meansub,'Color',plotcolor{i},'LineWidth',1);

%make significant section bold:
if sum(sigmask) ~= 0
    boldvals = meansub; 
    boldvals(~sigmask) = NaN; 
    plot(plottime, boldvals, 'Color',plotcolor{i}, 'LineWidth', 3.5); 
end

sigtime = time_values(sig_within);
if ~isempty(sigtime)
    for a = 1:length(sig_within)
        % Draw a horizontal line across the full x-axis at each significant value
        line([sigtime(a)-0.005 sigtime(a)+0.005], [y_pos(i) y_pos(i)], ...
             'Color', plotcolor{i}, 'LineWidth', 3);
    end
end

if ~isempty(sigtime)
%print first sig. time
x_pos = sigtime(1);
firstsigtime_label = round(sigtime(1)*1000,0); 
text(x_pos, y_pos(i)+0.005, [num2str(firstsigtime_label) 'ms'], 'FontSize', 10, 'Color',plotcolor{i});
end
grid 'on';
xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
ylabel('Balanced accuracy (chance=.5)');
% ylabel('Signed Dist. (a.u.)');
xlabel('time (s)');

xlim([-0.4 0.8]);
% ylim([0.45 0.6]);
title('Correct Superordinate Trials')

end


%% Stats balanced accuracies by PAS

y_pos(1) = 0.455;
y_pos(2) = 0.465;
y_pos(3) = 0.475;
y_pos(4) = 0.485;

plottitle{1} = 'Superordinate (N = 30)';
plottitle{2} = 'Basic (N = 30)';

timewindow = linspace(-0.4,0.8,241); %-0.4:0.01:0.8;
figure; set(gcf,'Color','w');
hold on;

for a = 1:2

subplot(1,2,a); hold on;
for pasrating = 1:4
pp = pasrating;
% if pasrating == 1
%     pp = 1;
% else
%     pp = 2:4;
% end
clearvars actual_dataset ds_group
for p=1:length(subjectslist)
    actual_dataset{p,1}.samples = nanmean(acc_for_stats{a}{p}(pp,:),1);
    actual_dataset{p,1}.sa.targets = 1;
    actual_dataset{p,1}.sa.chunks = p;
    actual_dataset{p,1}.fa.time = (1:241);
    actual_dataset{p,1}.fa.center_ids = (1:241);
    actual_dataset{p,1}.a.meeg.samples_field = 'trial';
    actual_dataset{p,1}.a.fdim.labels = {'time'};
    actual_dataset{p,1}.a.fdim.values{1,1} = linspace(-0.4,0.8,241);
end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group = cosmo_stack(ds_group);

allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group, 'time',allow_clustering_over_time);
ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);
sig_within = find(abs(ds_z.samples)>1.96);
sigmask = abs(ds_z.samples)>1.96;

clearvars temp
plottime = timewindow;

for s = 1:length(subjectslist)
    temp(s,:) = actual_dataset{s,1}.samples;
end

meansub = nanmean(temp);
sem_sub = nanstd(temp)/sqrt(length(subjectslist));
x = timewindow;
curve1 = meansub + sem_sub;
curve2 = meansub - sem_sub;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

plotcolor{1} = [7 26 187]./255;
plotcolor{2} = [44 187 164]./255; 
plotcolor{3} = [58 6 3]./255; 
plotcolor{4} = [234 54 128]./255;  

h = fill(x2, inBetween, plotcolor{pasrating});
set(h, 'EdgeColor', 'none'); % Remove black border
alpha(0.1);

time_values = timewindow;
plot1(pasrating) = plot(plottime,meansub,'Color',plotcolor{pasrating},'LineWidth',1);

%make significant section bold:
if sum(sigmask) ~= 0
    boldvals = meansub; 
    boldvals(~sigmask) = NaN; 
    plot(plottime, boldvals, 'Color',plotcolor{pasrating}, 'LineWidth', 3.5); 
end

sigtime = time_values(sig_within);
if ~isempty(sigtime)
    for i = 1:length(sig_within)
        % Draw a horizontal line across the full x-axis at each significant value
        line([sigtime(i)-0.005 sigtime(i)+0.005], [y_pos(pasrating) y_pos(pasrating)], ...
             'Color', plotcolor{pasrating}, 'LineWidth', 3);
    end
end

if ~isempty(sigtime)
%print first sig. time
x_pos = sigtime(1);
firstsigtime_label = round(sigtime(1)*1000,0); 
text(x_pos, y_pos(pasrating)+0.005, [num2str(firstsigtime_label) 'ms'], 'FontSize', 10, 'Color',plotcolor{pasrating});
end
grid 'on';
xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
ylabel('Balanced accuracy (chance=.5)');
xlabel('time (s)');

xlim([-0.4 0.8]);
% ylim([0.45 0.6]);
title(plottitle{a})
end
% legend([plot1(1) plot1(2) plot1(3) plot1(4)],{'PAS = 0','PAS = 1','PAS = 2','PAS = 3'});
end

%% Plot distances
timevals = linspace(-0.4,0.8,241);
smoothingval = 10;

summaryT = lmm_super_crosstime;

summaryT.DistZ = zeros(height(summaryT),1);
subjects = unique(summaryT.Subject);
traces   = unique(summaryT.SuperCorr);
for s = subjects'
    for tr = traces'
        idx = summaryT.Subject == s & summaryT.SuperCorr == tr;
        baseIdx = idx & summaryT.timebin <= 81;
        x = summaryT.Dist(idx);
        % mu = mean(x);
        % sigma = std(x);
        mu = nanmean(summaryT.Dist(baseIdx));
        sigma = nanstd(summaryT.Dist(baseIdx));

        summaryT.DistZ(idx) = (x - mu) / sigma;
    end
end

figure; set(gcf,'Color','w');
hold on;
subplot(1,2,1); hold on

colors = [
    94/255, 201/255,  98/255;
    68/255,   1/255,  84/255
];

groups = unique(summaryT.SuperCorr);
conds = unique(summaryT.Cond);
timeind = unique(summaryT.timebin);

for i = 1:length(groups)
    gval = groups(i);

    mean_vals = nan(size(timeind));
    sem_vals  = nan(size(timeind));

    for t = 1:length(timeind)
        idx = summaryT.SuperCorr == gval & summaryT.timebin == timeind(t);
        vals = summaryT.DistZ(idx);
        
        mean_vals(t) = mean(vals);
        sem_vals(t)  = std(vals) / sqrt(sum(~isnan(vals)));
    end

    shadedErrorBar(timevals, smoothdata(mean_vals,'sgolay',smoothingval), smoothdata(sem_vals,'sgolay',smoothingval), ...
        'lineProps', {'LineWidth', 2, 'DisplayName', sprintf('PAS = %d', gval), ...
                      'Color', colors(i,:)});
end

grid 'on';
xline(0,'--');
yline(0,'--');

xlabel('Time (s)');
ylabel('Signed Distance (a.u.)');
legend('Incorrect','Correct');
title('Superordinate distances');
xlim([-0.4,0.8]);
% ylim([-0.1 0.25]);

summaryT = lmm_basic_crosstime;

% subjects = unique(summaryT.Subject);
% Var_z = zeros(height(summaryT),1);
% for s = 1:length(subjects)
%     idx = (summaryT.Subject == subjects(s)) & (double(summaryT.Exemp) > 12);
%     v = summaryT.Dist(idx);
%     Var_z(idx) = (v - mean(v)) ./ std(v);
% 
%     idx2 = (summaryT.Subject == subjects(s)) & (double(summaryT.Exemp) <= 12);
%     v = summaryT.Dist(idx2);
%     Var_z(idx2) = (v - mean(v)) ./ std(v);
% end
% summaryT.Dist = Var_z;

summaryT.DistZ = zeros(height(summaryT),1);
subjects = unique(summaryT.Subject);
traces   = unique(summaryT.BasicCorr);
for s = subjects'
    for tr = traces'
        idx1 = summaryT.Subject == s & summaryT.BasicCorr == tr & ismember(double(summaryT.Exemp), 1:6);
        baseIdx = idx1 & summaryT.timebin <= 81;
        x = summaryT.Dist(idx1);
        % mu = mean(x);
        % sigma = std(x);
        mu = mean(summaryT.Dist(baseIdx));
        sigma = std(summaryT.Dist(baseIdx));

        summaryT.DistZ(idx1) = (x - mu) / sigma;
    end

    for tr = traces'
        idx2 = summaryT.Subject == s & summaryT.BasicCorr == tr & ismember(double(summaryT.Exemp), 7:12);
        baseIdx = idx2 & summaryT.timebin <= 81;
        x = summaryT.Dist(idx2);
        % mu = mean(x);
        % sigma = std(x);
        mu = nanmean(summaryT.Dist(baseIdx));
        sigma = nanstd(summaryT.Dist(baseIdx));

        summaryT.DistZ(idx2) = (x - mu) / sigma;
    end

end

subplot(1,2,2); hold on

groups = unique(summaryT.BasicCorr);
timeind = unique(summaryT.timebin);

for i = 1:length(groups)
    gval = groups(i);

    mean_vals = nan(size(timeind));
    sem_vals  = nan(size(timeind));

    for t = 1:length(timeind)
        idx = summaryT.BasicCorr == gval & summaryT.timebin == timeind(t);
        vals = summaryT.DistZ(idx);
        
        mean_vals(t) = nanmean(vals);
        sem_vals(t)  = nanstd(vals) / sqrt(sum(~isnan(vals)));
    end

    h = shadedErrorBar(timevals, smoothdata(mean_vals,'sgolay',smoothingval), smoothdata(sem_vals,'sgolay',smoothingval), ...
        'lineProps', {'LineWidth', 2, 'DisplayName', sprintf('PAS = %d', gval), ...
                      'Color', colors(i,:)});
    
    h.mainLine.DisplayName = sprintf('PAS = %d', gval); % for legend
end

grid 'on';
% ylim([-0.1 0.25]);
xline(0,'--');
yline(0,'--');
xlabel('Time (s)');
ylabel('Signed Distance (a.u.)');
legend('Correct','Incorrect');
title('Basic distances');
xlim([-0.4,0.8]);

%
