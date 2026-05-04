%% Bootstrapping onsets & peaks 10,000 times.

% This script takes the outputs from
% distancetobound_contentdecoding_bytrial.m and
% bootstrapps the decoding onset times.
% The bootstrapping procedure is the same as the other analysis from the
% MS: sample N subjects with replacement, repeat statistical test, and
% record the onset & peak times. Repeat 10,000 times to get a distribution
% for confidence intervals.

% Last updated July 1 2025, Ayaka Hachisuka

clear;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
%% Load each subject:
numSubjects = 33;
subjectslist = setdiff(1:numSubjects,[22 32 33]);
bootstrapNUM = 10000;
conditionLIST = {'AllCond'};
exemppairlist{1} = 1:3;
exemppairlist{2} = 4:6;
exemppairlist{3} = 7:9;
exemppairlist{4} = 10:12;

numCores = feature('numcores')-2;

for cc = 1:length(conditionLIST)
condition = conditionLIST{cc};
superorbasicLIST = {'Super','Basic'};
analysistypeLIST = {'Super','Basic'};

trialtype = 'ALL';
time_radius = 0;        
trial_bin_num = 1;
numrepeats = 100;
OCCIP = 0;
run_pca = 0;

for aa = 1:length(superorbasicLIST)
    superorbasic = superorbasicLIST{aa};
    analysistype = analysistypeLIST{aa};
    clearvars group_mvpaoutput
    sub_counter = 0;
for s = subjectslist
    sub_counter = sub_counter + 1;
    subjNum = s;
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end
    fprintf('Loading subject number: %d for %s category in %s\n',s,superorbasic,condition);

    if contains(superorbasic,'Super')
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/' ...
        'DistancetoBound_distbytrial_eucdistv4_withrepeat/Super/',condition,'/noPCA_35Hz_', num2str(time_radius), ...
        'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
         group_mvpaoutput{sub_counter} = load([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput'],'sl_map2');
    elseif contains(superorbasic,'Basic')
        save_dir1=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/' ...
    'DistancetoBound_distbytrial_eucdistv4_withrepeat/Basic1/',condition,'/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
        save_dir2=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/' ...
    'DistancetoBound_distbytrial_eucdistv4_withrepeat/Basic2/',condition,'/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_LDAba_',trialtype,'Trials/']);
        group_mvpaoutput{1,s} = load([save_dir1,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput'],'sl_map2');
        group_mvpaoutput{2,s} = load([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput'],'sl_map2');
    end
end

if size(group_mvpaoutput,1) == 2
    A = group_mvpaoutput{1,31}.sl_map2.dist;
    [max_val, linear_idx] = max(abs(A(:)));
    [row_idx, col_idx] = ind2sub(size(A), linear_idx);

    group_mvpaoutput{1,31}.sl_map2.dist(row_idx,:) = [];
    group_mvpaoutput{1,31}.sl_map2.samples(row_idx,:) = [];
    group_mvpaoutput{1,31}.sl_map2.origtrialinfo(row_idx,:) = [];
end

counter = 0;
for s = 1:size(group_mvpaoutput,2)
    %remove empties
    subjgrp = group_mvpaoutput(:,s);
    subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
    if ~isempty(subjgrp)
        counter = counter + 1;
        submatrix = cell2mat(cellfun(@(x) x.sl_map2.samples, subjgrp', 'UniformOutput', false)');
        avg{1,counter}.samples = nanmean(submatrix,1);
    end
end
n_subj = length(avg);
group_mvpaoutput2 = avg;

%% Group datasets together:
%some "non-unique" elements due to 0 being jittered, just use the
%parameters from subject 1 because they are mostly the same.
clearvars ds_group actual_dataset
actual_dataset = cell(n_subj,1);
for p=1:n_subj
    actual_dataset{p,1}.samples  = nanmean(group_mvpaoutput2{p}.samples,1);
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

%% Bootstrapping

timewindow = -0.4:0.005:0.8;
clearvars temp
plottime = timewindow;

%bootstrapping for statistics:
allow_clustering_over_time = true;

bootstrap_sig_onset = zeros(1,bootstrapNUM);
clearvars bootstrapgroup
parfor (b = 1:bootstrapNUM, numCores)
% for b = 1:bootstrapNUM
    bootstrap_indices = randi(n_subj, [1, n_subj]); %sample from all subjects, with replacement
    bootstrapacc = ds_group.samples(bootstrap_indices,:);

    bootstrapgroup = [];
    bootstrapgroup.samples = bootstrapacc;
    bootstrapgroup.sa = ds_group.sa;
    bootstrapgroup.fa = ds_group.fa;
    bootstrapgroup.a = ds_group.a;
    bootstrapgroup.a.fdim.values{1,1} = ds_group.a.fdim.values{1,1};
    bootstrapgroup.fa.center_ids = ds_group.fa.center_ids;
    bootstrapgroup.fa.time = ds_group.fa.time;
    
    nbrhood = cosmo_cluster_neighborhood(bootstrapgroup, 'time',allow_clustering_over_time);
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

savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/Dist2Bound_LMMstats/BootstrappingOnset_',num2str(bootstrapNUM),'_LDA_',trialtype,'Trials_stats/',condition,'/'];
if ~exist(savedir,'file'); mkdir(savedir); end
save(fullfile([savedir,'bootstrap_sigval_',analysistype]),'bootstrapbycond');
save(fullfile([savedir,'bootstrap_balacc_data_',analysistype]),'bootstrap_data_all');
save(fullfile([savedir,'ds_z_',analysistype]),'ds_z');

end
end

%% Bootstrapped onsets:
% Data saved from C2F_bootstrap_peakonset_groupfinal_supervsbasic.m
bootstrapNUM = 10000;
savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/Dist2Bound_LMMstats/BootstrappingOnset_',num2str(bootstrapNUM),'_LDA_',trialtype,'Trials_stats/',condition,'/'];

condColor{1} = [50 130 246]./255;
condColor{2} = [240 134 80]./255;

analysistypeLIST = {'Super','Basic'};
timewindow = linspace(-0.4,0.8,241);

clearvars bootstrapbycond_all
counter = 0;
for c = 1
    condition = conditionLIST{c};
for aa = 1:length(analysistypeLIST)
    counter = counter + 1;
    analysistype = analysistypeLIST{aa};
    temp = load(fullfile(savedir,['/ds_z_',analysistype,'.mat']),'ds_z');
    ds_z_all = temp.ds_z;
    for i = 1:length(ds_z_all)
        sigvals = find(abs(ds_z_all{i}.samples)>1.96);
        if ~isempty(sigvals)
            sig_onset(i) = timewindow(sigvals(1));
        else
            sig_onset(i) = nan;
        end
        [~,maxval] = max(ds_z_all{i}.samples);
        peakonset(i) = timewindow(maxval);
    end
    sig_onset_all{c,aa} = sig_onset;
    peak_all{c,aa} = peakonset;
    onsettime(c,aa) = nanmean(sig_onset);
    peaktime(c,aa) = nanmean(peakonset);
    onset_confidence_int(aa,:) = prctile(sig_onset, [2.5 97.5]);
    peak_confidence_int(aa,:) = prctile(peakonset, [2.5 97.5]); 

end
CIbyCond_onset{c} = onset_confidence_int; 
CIbyCond_peak{c} = peak_confidence_int; 

end

figure; set(gcf,'Color','w'); hold on
subplot(1,2,1); hold on
bar(1,onsettime(c,1),'FaceColor',condColor{1});
bar(2,onsettime(c,2),'FaceColor',condColor{2})
errorbar(1,onsettime(c,1), onsettime(c,1) - CIbyCond_onset{c}(1,1),CIbyCond_onset{c}(1,2) - onsettime(c,1),'.','Color','k','Linewidth',1);
errorbar(2,onsettime(c,2), onsettime(c,2) - CIbyCond_onset{c}(2,1),CIbyCond_onset{c}(2,2) - onsettime(c,2),'.','Color','k','Linewidth',1);
ylim([0 0.3])
xticks(1:2);
xticklabels({'Super','Basic'});
xtickangle(90);
ylabel('Onset Latency (ms)');
title(conditionLIST{c});
% set(gca, 'XDir', 'reverse'); % Flips the x-axis

subplot(1,2,2); hold on
hold on
bar(1,peaktime(c,1),'FaceColor',condColor{1});
bar(2,peaktime(c,2),'FaceColor',condColor{2})
errorbar(1,peaktime(c,1), peaktime(c,1) - CIbyCond_peak{c}(1,1),CIbyCond_peak{c}(1,2) - peaktime(c,1),'.','Color','k','Linewidth',1);
errorbar(2,peaktime(c,2), peaktime(c,2) - CIbyCond_peak{c}(2,1),CIbyCond_peak{c}(2,2) - peaktime(c,2),'.','Color','k','Linewidth',1);
ylim([0 0.3])
xticks(1:2);
xticklabels({'Super','Basic'});
xtickangle(90);
ylabel('Peak Latency (ms)');
title('Peak CI');

for c = 1:length(conditionLIST)
    disp(conditionLIST{c});
    disp('Onset - Super');
    fprintf('%.4f [%.4f %.4f]\n',  onsettime(c,1), CIbyCond_onset{c}(1,1), CIbyCond_onset{c}(1,2));
    disp('Onset - Basic');
    fprintf('%.4f [%.4f %.4f]\n',  onsettime(c,2), CIbyCond_onset{c}(2,1), CIbyCond_onset{c}(2,2));
    disp('Peak - Super');
    fprintf('%.4f [%.4f %.4f]\n',  peaktime(c,1), CIbyCond_peak{c}(1,1), CIbyCond_peak{c}(1,2));
    disp('Peak - Basic');
    fprintf('%.4f [%.4f %.4f]\n',  peaktime(c,2), CIbyCond_peak{c}(2,1), CIbyCond_peak{c}(2,2));
end
