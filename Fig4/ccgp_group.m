%% CCGP Group Analysis

% Train on Dog vs. Car, Test on Cat vs. Truck.
% Train on Dog vs. Truck, Test on Cat vs. Car.
% Train on Cat vs. Truck, Test on Dog vs. Car.
% Train on Cat vs. Car, Test on Dog vs. Truck.

%Last updated June 20 2025, Ayaka Hachisuka

clear;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
warning('off');
conditionLIST = {'Normal','Masked','LSF'};
trialtype = 'ALL';
num_subjects = 31;
n_modalities = 2; %it's a 2x2 cross-decoding schematic.
numrepeats =100;
timewindow = linspace(-0.4,0.8,241);

figure; set(gcf,'Color','w');
hold on
for c = 1:length(conditionLIST)
    condition = conditionLIST{c};
    savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/resultswith200Hzdata/',condition,'_BasicABCD', trialtype]);

counter1 = 0; counter2 = 0;
subjectslist = setdiff(1:num_subjects,[22 32 33]);

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

%% Average within-decoding subjects:
counter = 0;

for s = 1:size(sl_result_group_within,1)
    %remove empties
    subjgrp = sl_result_group_within(s,:);
    subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
    if ~isempty(subjgrp)
        counter = counter + 1;
        submatrix = cell2mat(cellfun(@(x) x.ds_searchlight_result.samples, subjgrp, 'UniformOutput', false)');
        avg_within{1,counter}.samples = mean(submatrix,1);
        avg_within{1,counter}.a.meeg.samples_field = 'trial';
        avg_within{1,counter}.a.fdim.labels = {'time'};
        avg_within{1,counter}.a.fdim.values{1,1} = linspace(-0.4,0.8,241);
        avg_within{1,counter}.fa.time = (1:241);
        avg_within{1,counter}.fa.center_ids = (1:241);
        avg_within{1,counter}.sa.targets = 1;
        avg_within{1,counter}.sa.chunks = s;
    end
end

%% Average across-decoding subjects:
counter = 0;
for s = 1:size(sl_result_group_across,1)
    %remove empties
    subjgrp = sl_result_group_across(s,:);
    subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
    if ~isempty(subjgrp)
        counter = counter + 1;
        submatrix = cell2mat(cellfun(@(x) x.ds_searchlight_result.samples, subjgrp, 'UniformOutput', false)');
        avg_across{1,counter}.samples = mean(submatrix,1);
        avg_across{1,counter}.a.meeg.samples_field = 'trial';
        avg_across{1,counter}.a.fdim.labels = {'time'};
        avg_across{1,counter}.a.fdim.values{1,1} = linspace(-0.4,0.8,241);
        avg_across{1,counter}.fa.time = (1:241);
        avg_across{1,counter}.fa.center_ids = (1:241);
        avg_across{1,counter}.sa.targets = 1;
        avg_across{1,counter}.sa.chunks = s;
    end
end
%% Group datasets together:
%some "non-unique" elements due to 0 being jittered, just use the
%parameters from subject 1 because they are mostly the same.

%within
for p=1:size(avg_within,2)
        actual_dataset{p,1}  = avg_within{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        actual_dataset{p,1}.a.fdim.values{1,1} = round(avg_within{p}.a.fdim.values{1,1},2);
end
[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group_within = cosmo_stack(ds_group);
[~,maxind] = max(nanmean(ds_group_within.samples));
realpeaktime1(c) = timewindow(maxind);

allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group_within, 'time',allow_clustering_over_time,'chan',false,'t',1);
ds_z = cosmo_montecarlo_cluster_stat(ds_group_within,nbrhood,'h0_mean',0.5,'niter',10000,'cluster_stat','maxsize','p_uncorrected',0.05);

sig_within = find(abs(ds_z.samples)>1.96);
sig_within_mask = abs(ds_z.samples)>1.96;
%across
for p=1:size(avg_across,2)
        actual_dataset{p,1}  = avg_across{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        actual_dataset{p,1}.a.fdim.values{1,1} = round(avg_across{p}.a.fdim.values{1,1},2);
end
[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group_across = cosmo_stack(ds_group);
[~,maxind] = max(nanmean(ds_group_across.samples));
realpeaktime2(c) = timewindow(maxind);

allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group_across, 'time',allow_clustering_over_time,'chan',false,'t',1);
ds_z = cosmo_montecarlo_cluster_stat(ds_group_across,nbrhood,'h0_mean',0.5,'niter',10000,'cluster_stat','maxsize','p_uncorrected',0.05);
sig_across = find(abs(ds_z.samples)>1.96);
sig_across_mask = abs(ds_z.samples)>1.96;
allbalacc{c,1} = ds_group_within.samples;
allbalacc{c,2} = ds_group_across.samples;
%% PLOTTING

condColor{1} = [155 155 155]./255;
condColor{2} = [0 18 154]./255;

clearvars temp
plottime = timewindow;
indstart = find(round(avg_within{1,1}.a.fdim.values{1,1},2) == plottime(1));
indend = find(round(avg_within{1,1}.a.fdim.values{1,1},2) == plottime(end));
inststart = indstart(1);
indend = indend(end);
%%% within
for s = 1:length(avg_within)
    temp(s,:) = avg_within{s}.samples(indstart:indend);
end
data_for_bf_within{c} = temp;

meansub = mean(temp);
sem_sub = std(temp)/sqrt(num_subjects);
x = timewindow;
curve1 = meansub + sem_sub;
curve2 = meansub - sem_sub;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

subplot(3,1,c);
hold on;

h = fill(x2, inBetween, condColor{1});
set(h, 'EdgeColor', 'none'); % Remove black border
alpha(0.3);

time_values = timewindow;

plot(plottime,meansub,'Color',condColor{1},'LineWidth',1);

%make significant section bold:
if sum(sig_within_mask) ~= 0; plot(plottime(sig_within_mask), meansub(sig_within_mask), 'Color',condColor{1}, 'LineWidth', 3.5); end

sigtime = time_values(sig_within);
for i = 1:length(sig_within)
    % Draw a horizontal line across the full x-axis at each significant value
    line([sigtime(i)-0.005 sigtime(i)+0.005], [0.48 0.48], ...
         'Color', condColor{1}, 'LineWidth', 3);
end

xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
ylabel('Balanced accuracy (chance=.5)');
% xlabel('time (s)');
% title(sprintf([condition, ' - CCGP, ' trialtype, ' (N = ', num2str(length(avg_within)), ')']));
title(sprintf([condition, ' (N = ', num2str(length(avg_within)), ')']));
xlim([-0.4 0.8]);
ylim([0.475 0.56]);
grid 'on';
%print first sig. time
if ~isempty(sigtime)
x_pos = sigtime(1);
y_pos = max(meansub)+0.03;
firstsigtime_label(c,1) = round(sigtime(1)*1000,0); 
text(x_pos, 0.485, [num2str(firstsigtime_label(c,1)) 'ms'], 'FontSize', 10, 'Color', condColor{1});
% ylim([0.48 0.55]);
end

%%% across
for s = 1:length(avg_across)
    temp(s,:) = avg_across{s}.samples(indstart:indend);
end
data_for_bf_across{c} = temp;

meansub = mean(temp);
sem_sub = std(temp)/sqrt(num_subjects);
x = timewindow;
curve1 = meansub + sem_sub;
curve2 = meansub - sem_sub;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

h = fill(x2, inBetween, condColor{2});
set(h, 'EdgeColor', 'none'); % Remove black border
alpha(0.3);

time_values = timewindow;

plot(plottime,meansub,'Color',condColor{2},'LineWidth',1);
grid 'on';
%make significant section bold:
if sum(sig_across_mask) ~= 0; plot(plottime(sig_across_mask), meansub(sig_across_mask), 'Color',condColor{2}, 'LineWidth', 3.5); end

sigtime = time_values(sig_across);

if ~isempty(sigtime)
    sigtime = time_values(sig_across);
    for i = 1:length(sig_across)
        % Draw a horizontal line across the full x-axis at each significant value
        line([sigtime(i)-0.005 sigtime(i)+0.005], [0.49 0.49], ...
             'Color', condColor{2}, 'LineWidth', 3);
    end
    
    %print first sig. time
    x_pos = sigtime(1);
    y_pos = max(meansub)+0.04;
    firstsigtime_label(c,2) = round(sigtime(1)*1000,0); 
    text(x_pos, 0.495, [num2str(firstsigtime_label(c,2)) 'ms'], 'FontSize', 10, 'Color', condColor{2});
    
end

%% cluster stats for differences
diff_acc = allbalacc{c,1} - allbalacc{c,2};
ds_group_within.samples = diff_acc;
allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group_within, 'time',allow_clustering_over_time);
ds_z = cosmo_montecarlo_cluster_stat(ds_group_within,nbrhood,'h0_mean',0,'niter',1000);
sig_ind = find(abs(ds_z.samples)>1.96);

corrsigtimeall = plottime(sig_ind);
if ~isempty(corrsigtimeall)
for i = 1:length(corrsigtimeall)
line([corrsigtimeall(i)-0.005 corrsigtimeall(i)+0.005], [0.55 0.55], ...
         'Color', 'k', 'LineWidth', 3);
end

diffsig_label = round(corrsigtimeall(1)*1000,0); 
x_pos = plottime(sig_ind(1));
text(x_pos, 0.544, [num2str(diffsig_label) 'ms'], 'FontSize', 10);

% legend({'','Within','','','Across'});
set(gca, 'FontName', 'Arial', 'FontSize', 18)
% legend boxoff  
end

%% BF testing 
% addpath(genpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/BFF_repo-master/'));
% [bf_within{c},bf_within_complement{c}] = bayesfactor(data_for_bf_within','nullvalue',0.5);
% [bf_across{c},bf_across_complement{c}] = bayesfactor(data_for_bf_across','nullvalue',0.5);

end

%% Bootstrapped onsets:
% Data saved from bootstrap_peakonset_groupfinal_supervsbasic.m

bootstrapNUM = 20000;
analysistypeLIST = {'Within','Across'};
timewindow = linspace(-0.4,0.8,241);

clearvars bootstrapbycond_all
counter = 0;
for c = 1:3
    condition = conditionLIST{c};
for aa = 1:length(analysistypeLIST)
    counter = counter + 1;
    analysistype = analysistypeLIST{aa};
    temp = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/',analysistype,'/BootstrappingTwoTail',num2str(bootstrapNUM),'_LDA_',trialtype,'Trials_stats/',condition, '/ds_z_',trialtype]),'ds_z');
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

% compute bootstrapped differences
onsetdiff = sig_onset_all{c,1} - sig_onset_all{c,2};
peakdiff = peak_all{c,1} - peak_all{c,2};

onsetdiff_ci{c} =  prctile(onsetdiff, [2.5 97.5]);
peakdiff_ci{c} =  prctile(peakdiff, [2.5 97.5]);

onsetdiff(isnan(onsetdiff)) = 0;

end

figure; set(gcf,'Color','w'); hold on
for c = 1:3
subplot(1,3,c); hold on
bar(1,onsettime(c,1),'FaceColor',condColor{1});
bar(2,onsettime(c,2),'FaceColor',condColor{2})
errorbar(1,onsettime(c,1), onsettime(c,1) - CIbyCond_onset{c}(1,1),CIbyCond_onset{c}(1,2) - onsettime(c,1),'.','Color','k','Linewidth',1);
errorbar(2,onsettime(c,2), onsettime(c,2) - CIbyCond_onset{c}(2,1),CIbyCond_onset{c}(2,2) - onsettime(c,2),'.','Color','k','Linewidth',1);
xlim([0 5]);
ylim([-0.4 0.8])
xticks(1:3);
xticklabels({'Control CCGP','CCGP'});
xtickangle(90);
ylabel('Onset Latency (ms)');
title(conditionLIST{c});
% set(gca, 'XDir', 'reverse'); % Flips the x-axis
end

figure; set(gcf,'Color','w'); hold on
for c = 1:3
subplot(1,3,c); hold on
bar(1,peaktime(c,1),'FaceColor',condColor{1});
bar(2,peaktime(c,2),'FaceColor',condColor{2})
errorbar(1,peaktime(c,1), peaktime(c,1) - CIbyCond_peak{c}(1,1),CIbyCond_peak{c}(1,2) - peaktime(c,1),'.','Color','k','Linewidth',1);
errorbar(2,peaktime(c,2), peaktime(c,2) - CIbyCond_peak{c}(2,1),CIbyCond_peak{c}(2,2) - peaktime(c,2),'.','Color','k','Linewidth',1);
xlim([0 5]);
ylim([-0.4 0.8])
xticks(1:3);
xticklabels({'Control CCGP','CCGP'});
xtickangle(90);
ylabel('Peak Latency (ms)');
title(conditionLIST{c});
% set(gca, 'XDir', 'reverse'); % Flips the x-axis
end

for c = 1:length(conditionLIST)
    disp(conditionLIST{c});
    disp('Onset - within');
    fprintf('%.4f [%.4f %.4f]\n',  onsettime(c,1), CIbyCond_onset{c}(1,1), CIbyCond_onset{c}(1,2));
    disp('Onset - across');
    fprintf('%.4f [%.4f %.4f]\n',  onsettime(c,2), CIbyCond_onset{c}(2,1), CIbyCond_onset{c}(2,2));
    disp('Peak - within');
    fprintf('%.4f [%.4f %.4f]\n',  peaktime(c,1), CIbyCond_peak{c}(1,1), CIbyCond_peak{c}(1,2));
    disp('Peak - across');
    fprintf('%.4f [%.4f %.4f]\n',  peaktime(c,2), CIbyCond_peak{c}(2,1), CIbyCond_peak{c}(2,2));
end

%% Sanity check plots:

for c = 1:3
    condition = conditionLIST{c};
    counter = 0;
for aa = 1:length(analysistypeLIST)
    for aa2 = 1:length(analysistypeLIST)
        if aa ~= aa2
            counter = counter + 1;
            sig_onset_all{c,aa}(isnan(sig_onset_all{c,aa})) = 0;
            sig_onset_all{c,aa2}(isnan(sig_onset_all{c,aa2}))= 0;
            sig_onset_diff = sig_onset_all{c,aa} - sig_onset_all{c,aa2};
            x = sig_onset_diff;
            onsetdiff_all{c,counter} = x;
            onset_diff_ci{c,counter} = prctile(x,[2.5 97.5]);
        end
    end
end
end

figcounter = 0;
figure; set(gcf,'Color','w');
for c = 1:3
    for counter = 1:2
        figcounter = figcounter + 1;
        subplot(3,3,figcounter); hold on
        histogram(onsetdiff_all{c,counter});
        xline(onset_diff_ci{c,1}(1));
        xline(onset_diff_ci{c,1}(2));
    end
end

%% Load permuted bootstrap data:
% conditions are scrambled.
% clc
% close all

clearvars confidence_int
figure; set(gcf,'Color','w');
for c = 1:3
    condition = conditionLIST{c};
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/Null_Within_Across/BootstrappingTwoTail20000_LDA_',trialtype,'Trials_stats/','/',condition,'/bootstrap_diff_onset',trialtype]),'bootstrap_diff_onset');
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/Null_Within_Across/BootstrappingTwoTail20000_LDA_',trialtype,'Trials_stats/','/',condition,'/bootstrap_diff_peak',trialtype]),'bootstrap_diff_peak');

    bootstrap_diff_onset(isnan(bootstrap_diff_onset)) = 0; %set non-sig differences to 0
    bootstrap_diff_peak(isnan(bootstrap_diff_peak)) = 0; %set non-sig differences to 0

    confidence_int = prctile(bootstrap_diff_onset, [2.5 97.5]); %low
    CIbyCond_onset{c} = confidence_int;

    subplot(2,3,c)
    hold on
    histogram(bootstrap_diff_onset);
    realonsetdiff = firstsigtime_label(c,2)/1000 - firstsigtime_label(c,1)/1000;
    xline(realonsetdiff,'r');
    xline(confidence_int(1));
    xline(confidence_int(2));

    confidence_int = prctile(bootstrap_diff_peak, [2.5 97.5]); %low
    CIbyCond_peak{c} = confidence_int;

    subplot(2,3,c+3)
    hold on
    histogram(bootstrap_diff_peak);
    realpeakdiff = realpeaktime2(c) - realpeaktime1(c);
    xline(realpeakdiff,'r');
    xline(confidence_int(1));
    xline(confidence_int(2));

    pval_onset(c) = (sum(abs(bootstrap_diff_onset) >= abs(realonsetdiff)) + 1) / (length(bootstrap_diff_onset) + 1);
    pval_onset_onetail(c) = (sum(bootstrap_diff_onset >= realonsetdiff) + 1) / (length(bootstrap_diff_onset) + 1);
    pval_peak(c) = (sum(abs(bootstrap_diff_peak) >= abs(realpeakdiff)) + 1) / (length(bootstrap_diff_peak) + 1);
end
fprintf('Normal: Onset p-val = %.4f, Peak p-val = %.4f\n',  pval_onset(1), pval_peak(1));
fprintf('Masked: Onset p-val = %.4f, Peak p-val = %.4f\n',  pval_onset(2), pval_peak(2));
fprintf('LSF: Onset p-val = %.4f, Peak p-val = %.4f\n',  pval_onset(3), pval_peak(3));

%% Comparing confidence intervals to 0

% conditions are scrambled.
% clc
% close all
conditionList = {'Normal','Masked','LSF'};
superorbasicLIST1 = {'Super','Super','Basic'};
superorbasicLIST2 = {'Basic','Exemplar','Exemplar'};

clearvars confidence_int superorbasic
figure; set(gcf,'Color','w');

for c = 1:3
    clearvars bootstrap_diff_onset bootstrap_diff_peak
    condition = conditionList{c};
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/DiffConfInt_WithinAcross/BootstrappingTwoTail_noclusteracrosstime_20000_LDA_',trialtype,['Trials_stats_TFCE' ...
        '/'],condition,'/bootstrap_diff_onset',trialtype]),'bootstrap_diff_onset');
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/DiffConfInt_WithinAcross/BootstrappingTwoTail_noclusteracrosstime_20000_LDA_',trialtype,['Trials_stats_TFCE' ...
        '/'],condition,'/bootstrap_diff_peak',trialtype]),'bootstrap_diff_peak');

    bootstrap_diff_onset(isnan(bootstrap_diff_onset)) = 0; %set non-sig differences to 0
    bootstrap_diff_peak(isnan(bootstrap_diff_peak)) = 0; %set non-sig differences to 0

    confidence_int = prctile(bootstrap_diff_onset, [2.5 97.5]); %low
    CIbyCond_onset{c} = confidence_int;

    subplot(3,3,c)
    hold on
    histogram(bootstrap_diff_onset);
    
    xline(confidence_int(1));
    xline(confidence_int(2));

    confidence_int = prctile(bootstrap_diff_peak, [2.5 97.5]); %low
    CIbyCond_peak{c} = confidence_int;

    subplot(3,3,c+6)
    hold on
    histogram(bootstrap_diff_peak);

    xline(confidence_int(1));
    xline(confidence_int(2));
end

disp('Onset differences:')
disp(CIbyCond_onset);

disp('Peak differences:')
disp(CIbyCond_peak);

