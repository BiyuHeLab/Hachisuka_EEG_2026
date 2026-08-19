%% CCGP Searchlight group analysis

clear;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/');
%% Load each subject:
conditionLIST = {'Normal'};
time_chunks = [0.168 0.316];
% conditionLIST = {'Masked'};
analysistype = 'ALL';
crosstype = {'Within','Across'};
numrepeats = 100;
for c = 1:length(conditionLIST)
    condition = conditionLIST{c};
    paths_C2F_EEG = C2F_EEG_SetPaths;
    vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
    vars = vars_C2F_EEG;
    time_radius = 0;        
    trial_bin_num = 4;
    n_modalities = 2;
    savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/resultswith200Hzdata/Searchlight_4neighbors/',condition,'_BasicABCD', analysistype]);
    numSubjects = 31;
    subjectslist = setdiff(1:numSubjects,[22 32 33]);
    scounter = 0 ;
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

%% Group datasets together:
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
        avg_within{1,counter}.a = subjgrp{1,1}.ds_searchlight_result.a;
        avg_within{1,counter}.fa = subjgrp{1,1}.ds_searchlight_result.fa;
        avg_within{1,counter}.sa = subjgrp{1,1}.ds_searchlight_result.sa;
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
        avg_across{1,counter}.a = subjgrp{1,1}.ds_searchlight_result.a;
        avg_across{1,counter}.fa = subjgrp{1,1}.ds_searchlight_result.fa;
        avg_across{1,counter}.sa = subjgrp{1,1}.ds_searchlight_result.sa;
    end
end
%% Group datasets together:
%some "non-unique" elements due to 0 being jittered, just use the
%parameters from subject 1 because they are mostly the same.
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
%within
for p=1:size(avg_within,2)
        actual_dataset{p,1}  = avg_within{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        % actual_dataset{p,1}.a.fdim.values{1,1} = round(avg_within{p}.a.fdim.values{1,1},2);
end
[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group_within = cosmo_stack(ds_group);

%across
for p=1:size(avg_across,2)
        actual_dataset{p,1}  = avg_across{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        % actual_dataset{p,1}.a.fdim.values{1,1} = round(avg_across{p}.a.fdim.values{1,1},2);
end
[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group_across = cosmo_stack(ds_group);

ds_group_all{1} = ds_group_within;
ds_group_all{2} = ds_group_across;

%% prepare layout

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/forward');
data_tl.label = actual_dataset{1,1}.a.fdim.values{1,1};
easycaplayout = load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/layout_acticap128chan.mat').lay_acticap128chan;

common_labels = intersect(easycaplayout.label,data_tl.label);
[~, idx] = ismember(common_labels, easycaplayout.label);
new_layout = struct();

new_layout.pos = easycaplayout.pos(idx, :);
new_layout.width = easycaplayout.width(idx);
new_layout.height = easycaplayout.height(idx);
new_layout.label = easycaplayout.label(idx);
new_layout.outline = easycaplayout.outline;
new_layout.mask = easycaplayout.mask;
new_layout.cfg = easycaplayout.cfg;
layout = new_layout;
layout.cfg.channel = new_layout.label;

[~, sortIdx] = sort(layout.pos(:, 2)); % Sort by y-coordinate
sorted_electrodes = layout.label(sortIdx); % Sorted electrode names
sorted_positions = layout.pos(sortIdx, :); % Sorted positions
sorted_layout = layout;
sorted_layout.label = sorted_electrodes;
sorted_layout.pos = sorted_positions;

allow_clustering_over_time = true;

cfg.layout = new_layout;
cfg.method = 'triangulation';
lay = ft_prepare_neighbours(cfg)'; %cosmo_meeg_chan_neighbors should do the same as this funtion          

for i = 1:2
chan_nbrhood=cosmo_meeg_chan_neighborhood(ds_group_all{i}, lay);
nbrhood = cosmo_cluster_neighborhood(ds_group_all{i},'chan',chan_nbrhood,'time',allow_clustering_over_time);

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');
paths_C2F_EEG = C2F_EEG_SetPaths;
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
vars = vars_C2F_EEG;

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
% ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',10000);
ds_z = cosmo_montecarlo_cluster_stat(ds_group_all{i},nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);

sig = find(abs(ds_z.samples)>1.96);

%% Topoplots:

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/utilities');

ds_group_ft = cosmo_map2meeg(ds_group_all{i});
clearvars ds_z_ft
ds_z_ft = cosmo_map2meeg(ds_z);

avggroup = squeeze(nanmean(ds_group_ft.trial,1));

z_scores = ds_z_ft.avg;
xval = round(linspace(-0.4, 0.8, 241),2);
zeroval = find(xval == 0);

plotval = avggroup; %z_scores, avggroup
sortedplotval = plotval(sortIdx,:);
sorted_z = z_scores(sortIdx,:);
sortedlabel = ds_z_ft.label(sortIdx);
figure; set(gcf,'Color','w'); hold on;
imagesc(sortedplotval);
% Outline sig elements
contour(sorted_z, 1, 'LineColor', 'k', 'LineWidth', 1.5);

xticks(linspace(1, size(plotval, 2), 50)); % Choose the number of ticks you want
xticklabels(round(linspace(-0.4, 0.8, 50),2)); % Set custom labels from -0.5 to 0.5
xline(zeroval(1), '--', 'LineWidth', 1, 'Color', 'k');
xtickangle(45);
yticks(1:length(ds_z_ft.label));
yticklabels(sortedlabel);
ylabel('channel #')
xlabel('time(s)');
clim([0.48 0.55]);
set(gca, 'FontSize', 8); % Adjust the font size as needed
colorbar();
title(sprintf([condition,' - ' crosstype{i}]))

%%
channel_labels = ds_z_ft.label;  
time_vector = linspace(-0.4, 0.8, size(z_scores, 2)); 

% initalize variable:
significant_channels = cell(length(time_chunks) - 1, 1);

% Loop through each time chunk to find significant channels
for t = 1:length(time_chunks)
    time_idx = find(round(time_vector,3) == round(time_chunks(t),2));

    significant = any(abs(z_scores(:, time_idx)) > 1.96,2);  % Logical vector for significant channels
    significant_channels{t} = channel_labels(significant);
end

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/plotting');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/fileio');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/forward');
       
nColors = 256;
if contains(condition,'Normal')
    condColor = [0.9290, 0.6940, 0.1250];
elseif contains(condition,'Masked')
    condColor = [0.6350, 0.0780, 0.1840];
elseif contains(condition,'LSF')
    condColor = [0.4660 0.6740 0.1880];
end

gray_start = [1 1 1];

% Build colormap from gray to vivid color
custom_cmap = [linspace(gray_start(1), condColor(1), nColors)', ...
               linspace(gray_start(2), condColor(2), nColors)', ...
               linspace(gray_start(3), condColor(3), nColors)'];

figure; set(gcf,'Color','w');
hold on
for t = 1:length(time_chunks)
    ax = subplot(1, round(length(time_chunks)), t);  % Create a 3x3 grid of subplots for 9 chunks
    % cfg.xlim = [time_chunks(t), time_chunks(t + 1)];  % Time window for this subplot
    cfg.xlim = [time_chunks(t),time_chunks(t)];  % Time window for this subplot
    cfg.zlim=0.49:0.01:0.54;
    cfg.markerfontsize = 1;
    cfg.markercolor = [0.5 0.5 0.5];
    cfg.colormap = custom_cmap;
    cfg.style = 'both';
    cfg.highlight = 'on';
    cfg.highlightchannel = significant_channels{t};
    cfg.highlightsymbol = '.';
    cfg.highlightcolor = 'k';
    cfg.highlightsize = 12;
    cfg.comment = 'xlim';
    ft_topoplotER(cfg, ds_group_ft);
end
title(sprintf([condition,'-', crosstype{i}]))
set(gca, 'FontName', 'Arial', 'FontSize', 25)
end
end
