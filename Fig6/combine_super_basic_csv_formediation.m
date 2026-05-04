%% Plotting PAS vs. Distance to Bound (LDA)

% This script plots the following (with outputs from
% distancetobound_contentdecoding_bytrial.m)
% reformatted for mediation analysis.

% Last updated July 1 2025, Ayaka Hachisuka

clear; 
num_subjects = 31;
subjectslist = setdiff(1:num_subjects,[22 32 33]);

cond = 'AllCond';
trialtype = 'SuperCorrect';

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

%% Super
load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/Behavioral/trialinfo_alltrials_30subsV2.mat','trialinfo');
save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/' ...
    'DistancetoBound_distbytrial_eucdistv4_withrepeat/Super/',cond,'/noPCA_35Hz_', num2str(time_radius), ...
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

    sl_map = sl_map2;

    trialnums = sl_map.origtrialinfo(:,4);
    trialnums(isnan(trialnums)) = [];
    [sorted_trials, sortIdx] = sortrows(trialnums, 1);

    sorted_dist = sl_map.dist(sortIdx,:);

    sorted_pas1 = sl_map.origtrialinfo(:,19);
    sorted_pas = sorted_pas1(sortIdx,:);
    super_sortedpas = sorted_pas;

    sorted_exemptrials1 = sl_map.origtrialinfo(:,12);
    sorted_exemptrials = sorted_exemptrials1(sortIdx,:);

    sorted_cond = sl_map.origtrialinfo(:,9);
    sorted_cond = sorted_cond(sortIdx,:);

    trialnumbers = NaN(size(sorted_dist,1), 1);                
    valid_rows = any(~isnan(sorted_dist), 2); % rows without NaN
    trialnumbers(valid_rows) = find(valid_rows); % Add index is the row is valid.

    %remove trials where PAS is NaN (no response):
    pasresp = sorted_pas;
    condlist = sorted_cond;
    nantrials = find(isnan(pasresp));
    if ~isempty(nantrials)
        condlist(nantrials,:) = [];
        pasresp(nantrials,:) = [];
        sorted_dist(nantrials,:) = [];
        sorted_exemptrials(nantrials,:) = []; 
    end

    % SIGNED DISTANCE
    pos_basic_trials = 1:6;
    neg_basic_trials = 7:12;
    trialinfo2 = nan(size(sorted_exemptrials));
    trialinfo2(ismember(sorted_exemptrials,pos_basic_trials)) = 1; % dog trials
    trialinfo2(ismember(sorted_exemptrials, neg_basic_trials)) = 2; % cat trials
    trialinfo2_expanded = repmat(trialinfo2, 1, 241);
    % if it's an inanimate trial, multiply by -1. This keeps the negative
    % sign of misclassified distances.
    sorted_dist(trialinfo2_expanded == 2) = -1 * sorted_dist(trialinfo2_expanded == 2);

    % correct baseline
    baseline = nanmean(sorted_dist(:,1:81),2);
    sorted_dist = sorted_dist - baseline;

    alldist_lmm{s} = zscore(sorted_dist,0,1);
    allpas_lmm{s} = pasresp;
    allcond_lmm{s} = condlist;
    trialnum_lmm{s} = 1:length(pasresp);
end

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
        cond_lmm = [cond_lmm; allcond_lmm{s}];
        pas_lmm = [pas_lmm; allpas_lmm{s}];        
        dist_lmm = [dist_lmm; alldist_lmm{s}(:,t)];      
        trial_lmm = [trial_lmm; trialnum_lmm{s}'];       
        subj_lmm = [subj_lmm; repmat(s,length(allpas_lmm{s}),1)];       
    end
    pas_lmm(pas_lmm == 3) = 4;
    pas_lmm(pas_lmm == 2) = 3;
    pas_lmm(pas_lmm == 1) = 2;
    pas_lmm(pas_lmm == 0) = 1;

    lmm_super_table{t} = table(categorical(pas_lmm), dist_lmm, cond_lmm, trial_lmm, subj_lmm,'VariableNames', {'PAS', 'Dist','Cond','Trial', 'Subject'});
end
%% Basic 1
clearvars distbyexemp alldist_lmm
trialsofinterest = 1:18;
load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/Behavioral/trialinfo_alltrials_30subsV2.mat','trialinfo');
save_dir1=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/' ...
    'DistancetoBound_distbytrial_eucdistv4_withrepeat/Basic1/',cond,'/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
save_dir2=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/' ...
    'DistancetoBound_distbytrial_eucdistv4_withrepeat/Basic2/',cond,'/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
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
    basic1_sl_map = sl_map2;

    temp_origtrialinfo = basic1_sl_map.origtrialinfo;
    temp_origtrialinfo1 = temp_origtrialinfo;
    temp_origtrialinfo1(all(isnan(temp_origtrialinfo), 2), :) = [];
    basic1_sl_map.origtrialinfo = temp_origtrialinfo1;
    basic1_sl_map.dist(all(isnan(temp_origtrialinfo), 2), :) = [];

    if s == 30
        A = basic1_sl_map.dist;
        [max_val, linear_idx] = max(abs(A(:)));
        [row_idx, col_idx] = ind2sub(size(A), linear_idx);
    
        basic1_sl_map.dist(row_idx,:) = [];
        basic1_sl_map.samples(row_idx,:) = [];
        basic1_sl_map.origtrialinfo(row_idx,:) = [];
    end

    load([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',cond,'_mvpaoutput'],'sl_map2');
    basic2_sl_map = sl_map2;

    temp_origtrialinfo = basic2_sl_map.origtrialinfo;
    temp_origtrialinfo1 = temp_origtrialinfo;
    temp_origtrialinfo1(all(isnan(temp_origtrialinfo), 2), :) = [];
    basic2_sl_map.origtrialinfo = temp_origtrialinfo1;
    basic2_sl_map.dist(all(isnan(temp_origtrialinfo), 2), :) = [];

    combined_basic = [basic1_sl_map.origtrialinfo(:,4); basic2_sl_map.origtrialinfo(:,4)];
    combined_basic(isnan(combined_basic)) = [];
    [sorted_trials, sortIdx] = sortrows(combined_basic, 1);

    combined_basic_dist = [basic1_sl_map.dist; basic2_sl_map.dist];
    combined_basic_dist = combined_basic_dist(sortIdx,:);

    combined_basic_pas = [basic1_sl_map.origtrialinfo(:,19); basic2_sl_map.origtrialinfo(:,19)];
    combined_basic_pas = combined_basic_pas(sortIdx,:);
    basic_sortedpas = sorted_pas;

    combined_basic_exemptrials = [basic1_sl_map.origtrialinfo(:,12); basic2_sl_map.origtrialinfo(:,12)];
    combined_basic_exemptrials = combined_basic_exemptrials(sortIdx,:);

    combined_basic_cond = [basic1_sl_map.origtrialinfo(:,9); basic2_sl_map.origtrialinfo(:,9)];
    combined_basic_cond = combined_basic_cond(sortIdx,:);

    trialnumbers = NaN(size(combined_basic_dist,1), 1);                
    valid_rows = any(~isnan(combined_basic_dist), 2); % rows without NaN
    trialnumbers(valid_rows) = find(valid_rows); % Add index is the row is valid.

    %remove trials where PAS is NaN (no response):
    pasresp = combined_basic_pas;
    condlist = combined_basic_cond;
    nantrials = find(isnan(pasresp));
    if ~isempty(nantrials)
        pasresp(nantrials,:) = [];
        condlist(nantrials,:) = [];
        combined_basic_dist(nantrials,:) = [];
        combined_basic_exemptrials(nantrials,:) = []; 
    end

    % SIGNED DISTANCE
    pos_basic_trials = [1 2 3 4 7 8 9];
    neg_basic_trials = [4 5 6 10 11 12];
    trialinfo2 = nan(size(combined_basic_exemptrials));
    trialinfo2(ismember(combined_basic_exemptrials,pos_basic_trials)) = 1; % dog trials
    trialinfo2(ismember(combined_basic_exemptrials, neg_basic_trials)) = 2; % cat trials
    trialinfo2_expanded = repmat(trialinfo2, 1, 241);
    % if it's an inanimate trial, multiply by -1. This keeps the negative
    % sign of misclassified distances.
    combined_basic_dist(trialinfo2_expanded == 2) = -1 * combined_basic_dist(trialinfo2_expanded == 2);

    idx1 = trialinfo2 == 1;
    combined_basic_dist(idx1,:) = zscore(combined_basic_dist(idx1,:),0,1);

    idx2 = trialinfo2 == 2;
    combined_basic_dist(idx2,:) = zscore(combined_basic_dist(idx2,:),0,1);

    % correct baseline
    baseline = nanmean(combined_basic_dist(:,1:81),2);
    combined_basic_dist = combined_basic_dist - baseline;

    alldist_lmm{s} = combined_basic_dist;
    allpas_lmm{s} = pasresp;
    allcond_lmm{s} = condlist;
    trialnum_lmm{s} = 1:length(pasresp);
end


%%
lmm_basic = cell(241,1);
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
        cond_lmm = [cond_lmm; allcond_lmm{s}];
        pas_lmm = [pas_lmm; allpas_lmm{s}];        
        dist_lmm = [dist_lmm; alldist_lmm{s}(:,t)];      
        trial_lmm = [trial_lmm; trialnum_lmm{s}'];       
        subj_lmm = [subj_lmm; repmat(s,length(allpas_lmm{s}),1)];       
    end
    pas_lmm(pas_lmm == 3) = 4;
    pas_lmm(pas_lmm == 2) = 3;
    pas_lmm(pas_lmm == 1) = 2;
    pas_lmm(pas_lmm == 0) = 1;

    lmm_basic_table{t} = table(categorical(pas_lmm), dist_lmm, cond_lmm, trial_lmm, subj_lmm,'VariableNames', {'PAS', 'Dist','Cond','Trial', 'Subject'});
end
%% Re-organizing variables and saving them for CLMM:
for t = 1:241
    lmm_basic_table{t}.PAS = double(lmm_basic_table{t}.PAS);
end

lmm_super_crosstime = table();
for t = 1:241
    temp = lmm_super_table{t};
    temp.timebin = repmat(t, height(temp), 1);
    lmm_super_crosstime = [lmm_super_crosstime; temp];
end

lmm_basic_crosstime = table();
for t = 1:241
    temp = lmm_basic_table{t};
    temp.timebin = repmat(t, height(temp), 1);
    lmm_basic_crosstime = [lmm_basic_crosstime; temp];
end

% save data
writetable(lmm_super_crosstime, ['/isilon/LFMI/VMdrive/Ayaka/EEG/Distance2Bound_LMMstats/Signed_dist/',trialtype,'/lmm_superV3.csv']);
writetable(lmm_basic_crosstime, ['/isilon/LFMI/VMdrive/Ayaka/EEG/Distance2Bound_LMMstats/Signed_dist/',trialtype,'/lmm_basicV3.csv']);

%% 

figure; hold on
for p = 1:4
    subplot(2,2,p); hold on
for s = 1:30
    idx = (lmm_basic_crosstime.Subject == s) & (lmm_basic_crosstime.PAS == p);
    data = lmm_basic_crosstime(idx, :);
    bins = unique(data.timebin);
    meanTrace = zeros(length(bins), 1);
    for i = 1:length(bins)
        binIdx = data.timebin == bins(i);
        meanTrace(i) = mean(data.Dist(binIdx), 'omitnan');
    end
    if meanTrace ~= 0
    plot(bins, zscore(meanTrace), 'LineWidth', 1);
    meanPASdist{p}(s,:) = zscore(meanTrace);
    end
end
end

plottime = linspace(-0.4,0.8,241);
figure; hold on
for p = 1:4
    temp = mean(meanPASdist{p});
    plot(plottime,temp);
end
xline(0,'--');
xline(0.062,'--');
xline(0.112,'--');
yline(0,'--');
legend({'PAS=0','PAS=1','PAS=2','PAS=3'})
