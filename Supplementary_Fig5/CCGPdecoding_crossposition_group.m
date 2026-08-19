%% CCGP cross position decoding group

clear; close all
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/ADAM-1.14-beta/install');
% startup;
%% Load each subject:
numSubjects = 33;
subjectslist = setdiff(1:numSubjects,[22 32 33]);
numrepeats = 100;

conditionlist = {'Normal','Masked','LSF'};
% analysistype = 'CCGPCrossPos_leaveonetrialout';
% analysistype = 'CCGPByPos';
analysistype = 'CCGPCrossPos';
trialtype = 'ALL';


figure; set(gcf,'Color','w');
hold on

cc_counter = 0;

for c = 1:length(conditionlist)

condition = conditionlist{c};

clearvars avg_vals group_mvpaoutput sig sigtime
condition = conditionlist{c};
time_radius = 0;        
trial_bin_num = 4;
nchunks = 10;
OCCIP = 0;
run_pca = 0;
n_modalities = 2;

%plot colors
plotcolor= 'b';

counter1 = 0; counter2 = 0;

for posvalue1 = 1:4

% savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Resubmission_results/AllChan/',analysistype,'/',condition,'_BasicABCD', trialtype,'_pos',num2str(posvalue1)]);
savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Final_Resubmission_results_April16/Fig4/',analysistype,'/',condition,'_BasicABCD', trialtype,'_pos',num2str(posvalue1)]);

for s = subjectslist
    for pair = 1:2
        for train_modality=1:n_modalities
            for test_modality=1:n_modalities
                if train_modality == test_modality
                    filepathname = fullfile([savedir, '/sub', num2str(s), '_pair', num2str(pair), '_combo_', ...
                    num2str(train_modality), num2str(test_modality), '_',num2str(numrepeats),'repeats_sl_map.mat']);
                    if exist(filepathname,'file')
                        counter1 = counter1 + 1;
                        sl_result_group_within{s,counter1} = load(filepathname,'sl_map');
                    end
                else
                    filepathname = fullfile([savedir, '/sub', num2str(s), '_pair', num2str(pair), '_combo_', ...
                    num2str(train_modality), num2str(test_modality), '_',num2str(numrepeats),'repeats_sl_map.mat']);
                    if exist(filepathname,'file')
                        counter2 = counter2 + 1;
                        sl_result_group_across{s,counter2} = load(fullfile([savedir, '/sub', num2str(s), '_pair', num2str(pair), '_combo_', ...
                        num2str(train_modality), num2str(test_modality), '_',num2str(numrepeats),'repeats_sl_map']),'sl_map');
                    end
                end
            end
        end
    end
    end
end

counter = 0;
for s = 1:size(sl_result_group_within,1)
    %remove empties
    subjgrp = sl_result_group_within(s,:);
    subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
    if ~isempty(subjgrp)
        counter = counter + 1;
        submatrix = cell2mat(cellfun(@(x) x.sl_map.samples, subjgrp, 'UniformOutput', false)');
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
        submatrix = cell2mat(cellfun(@(x) x.sl_map.samples, subjgrp, 'UniformOutput', false)');
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

timewindow = linspace(-0.4,0.8,241);
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
% ds_z = cosmo_montecarlo_cluster_stat(ds_group_within,nbrhood,'h0_mean',0.5,'niter',1000);
ds_z = cosmo_montecarlo_cluster_stat(ds_group_within,nbrhood,'h0_mean',0.5,'niter',10000,'cluster_stat','maxsize','p_uncorrected',0.05);
% sig_within = find(ds_z.samples>1.6449);
sig_within = find(abs(ds_z.samples)>1.96);
sig_within_mask = abs(ds_z.samples)>1.96;

%% across
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
% ds_z = cosmo_montecarlo_cluster_stat(ds_group_across,nbrhood,'h0_mean',0.5,'niter',1000);
ds_z = cosmo_montecarlo_cluster_stat(ds_group_across,nbrhood,'h0_mean',0.5,'niter',10000,'cluster_stat','maxsize','p_uncorrected',0.05);
sig_across = find(abs(ds_z.samples)>1.96);
sig_across_mask = abs(ds_z.samples)>1.96;
allbalacc{c,1} = ds_group_within.samples;
allbalacc{c,2} = ds_group_across.samples;
%% PLOTTING

condColor{1} = [155 155 155]./255; %[82 202 255]./255;
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
sem_sub = std(temp)/sqrt(numSubjects);
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
[val, ind] = max(meansub);
subjlevel_peak(c,1) = time_values(ind);

%make significant section bold:
% if sum(sig_within_mask) ~= 0; plot(plottime(sig_within_mask), meansub(sig_within_mask), 'Color',condColor{1}, 'LineWidth', 3.5); end
%make significant section bold:

bold_meansub = meansub;
bold_meansub(~sig_within_mask) = NaN; % set mask to NaN where it is 0
%make significant section bold:
sigtime = time_values(sig_within);
if sum(sig_within_mask) ~= 0; plot(plottime, bold_meansub, 'Color',condColor{1}, 'LineWidth', 3.5); end
for i = 1:length(sigtime)
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
subjlevel_onset(c,1) = sigtime(1);
text(x_pos, 0.485, [num2str(firstsigtime_label(c,1)) 'ms'], 'FontSize', 10, 'Color', condColor{1});
% ylim([0.48 0.55]);
end

%%% across
for s = 1:length(avg_across)
    temp(s,:) = avg_across{s}.samples;
end
data_for_bf_across{c} = temp;

meansub = nanmean(temp);
sem_sub = nanstd(temp)/sqrt(numSubjects);
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
[val, ind] = max(meansub);
subjlevel_peak(c,2) = time_values(ind);

grid 'on';
%make significant section bold:
% if sum(sig_across_mask) ~= 0; plot(plottime(sig_across_mask), meansub(sig_across_mask), 'Color',condColor{2}, 'LineWidth', 3.5); end
bold_meansub = meansub;
bold_meansub(~sig_across_mask) = NaN; % set mask to NaN where it is 0
%make significant section bold:
sigtime = time_values(sig_across);
if sum(sig_across_mask) ~= 0; plot(plottime, bold_meansub, 'Color',condColor{2}, 'LineWidth', 3.5); end
for i = 1:length(sigtime)
    line([sigtime(i)-0.005 sigtime(i)+0.005], [0.48 0.48], ...
         'Color', condColor{2}, 'LineWidth', 3);
end

clearvars sigtime
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
    subjlevel_onset(c,2) = sigtime(1);
    text(x_pos, 0.495, [num2str(firstsigtime_label(c,2)) 'ms'], 'FontSize', 10, 'Color', condColor{2});
    
end

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
end

ax = gca;
ax.XTick = -0.4:0.2:0.8;
ax.XMinorTick = 'on';
ax.XAxis.MinorTickValues = -0.4:0.1:0.8;

% legend({'','Within','','','Across'});
set(gca, 'FontName', 'Arial', 'FontSize', 18)
% legend boxoff  

end
