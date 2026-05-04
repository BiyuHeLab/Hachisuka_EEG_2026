%% MEEG timelock searchlight group analysis
% Last updated Feb 6 2026, Ayaka Hachisuka

clear;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');
%% Load each subject:
conditionLIST = {'Normal','Masked','LSF'};
% conditionLIST = {'LSF'};
time_chunks = [0.05 0.06 0.07 0.08 0.09 0.1 0.3 0.35 0.36];
% time_chunks = [0 0.22 0.35];
% time_chunks = [0 0.13 0.2];

% time_chunks = [0 0.16 0.23];
% time_chunks = [0 0.12 0.2];

% time_chunks = [0 0.1 0.22];
% time_chunks = [0 0.11 0.21];
for c = 1:length(conditionLIST)
    condition = conditionLIST{c};
    superorbasic = 'Super';
    time_radius = 0;        
    trial_bin_num = 4;

    if strcmp(superorbasic,'Super')
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/Searchlight/noPCA_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan',superorbasic,'_LDA_AllTrials/']);
    elseif strcmp(superorbasic,'Basic')
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/Searchlight/noPCA_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan',superorbasic,'1_LDA_AllTrials/']);
        save_dir2=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/Searchlight/noPCA_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan',superorbasic,'2_LDA_AllTrials/']);
    end
    
    numSubjects = 31;
    subjectslist = setdiff(1:numSubjects,[22 32 33]);
    scounter = 0;
for s = subjectslist
    scounter = scounter + 1;
    subjNum = s;
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end
    
    mvpaoutputdir = fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'sl_tl_ds.mat']); % this is the direct MVPA output  
    group_mvpaoutput{scounter} = load(mvpaoutputdir).sl_tl_ds;
    %%%%%%%%%
    if exist('save_dir2','var')
        mvpaoutputdir = fullfile([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',condition,'sl_tl_ds.mat']); % this is the direct MVPA output
        group_mvpaoutput2{scounter} = load(mvpaoutputdir).sl_tl_ds;
    end
end
scounter = 0;
if exist('save_dir2','var')
    for s = subjectslist
        scounter = scounter + 1;
        group_mvpaoutput_avg{1,scounter}.samples = mean([group_mvpaoutput{1,scounter}.samples; group_mvpaoutput2{1,scounter}.samples]);
        group_mvpaoutput_avg{1,scounter}.a = group_mvpaoutput{1,scounter}.a;
        group_mvpaoutput_avg{1,scounter}.fa = group_mvpaoutput{1,scounter}.fa;
        group_mvpaoutput_avg{1,scounter}.sa = group_mvpaoutput{1,scounter}.sa;
    end
    group_mvpaoutput = group_mvpaoutput_avg;
end

n_subj = length(group_mvpaoutput);

%% Group datasets together:
%some "non-unique" elements due to 0 being jittered, just use the
%parameters from subject 1 because they are mostly the same.

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
n_subj = length(group_mvpaoutput);
ds_cell = cell(n_subj,1);
for p=1:n_subj
        actual_dataset{p,1}  = group_mvpaoutput{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        %round values to the 3rd decimal because "0" is a very small number
        %that jitters around 0.
        actual_dataset{p,1}.a.fdim.values{2,1} = round(group_mvpaoutput{p}.a.fdim.values{2,1},3);
end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group = cosmo_stack(ds_group);


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

chan_nbrhood=cosmo_meeg_chan_neighborhood(ds_group, lay);
nbrhood = cosmo_cluster_neighborhood(ds_group,'chan',chan_nbrhood,'time',allow_clustering_over_time);

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');
paths_C2F_EEG = C2F_EEG_SetPaths;
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
vars = vars_C2F_EEG;

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.01);

sig = find(abs(ds_z.samples)>1.96);

%% Topoplots:

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/utilities');

ds_group_ft = cosmo_map2meeg(ds_group);
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
title(sprintf([condition,'-',superorbasic]))

%%
channel_labels = ds_z_ft.label;  
time_vector = linspace(-0.4, 0.8, size(z_scores, 2)); 

% initalize variable:
significant_channels = cell(length(time_chunks) - 1, 1);

% Loop through each time chunk to find significant channels
for t = 1:length(time_chunks)
    % time_idx = find(time_vector >= time_chunks(t) & time_vector < time_chunks(t + 1));
    time_idx = find(round(time_vector,3) == time_chunks(t));

    significant = any(abs(z_scores(:, time_idx)) > 1.96,2);  % Logical vector for significant channels
    significant_channels{t} = channel_labels(significant);
end

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/plotting');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/fileio');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/forward');
cfg.layout = layout;   
nColors = 256;
if contains(condition,'Normal')
    condColor = [0.9290, 0.6940, 0.1250];
elseif contains(condition,'Masked')
    condColor = [0.6350, 0.0780, 0.1840];
elseif contains(condition,'LSF')
    condColor = [0.4660 0.6740 0.1880];
end
gray_start = [1 1 1];

% build colormap from gray to vivid color
custom_cmap = [linspace(gray_start(1), condColor(1), nColors)', ...
               linspace(gray_start(2), condColor(2), nColors)', ...
               linspace(gray_start(3), condColor(3), nColors)'];

addpath(pwd);
figure; set(gcf,'Color','w');
hold on
for t = 1:length(time_chunks)
    ax = subplot(1,length(time_chunks),t);
    cfg.xlim = [time_chunks(t),time_chunks(t)];  % Time window for this subplot
    cfg.zlim=0.49:0.01:0.54;
    cfg.markerfontsize = 1;
    cfg.markercolor = [0.5 0.5 0.5];
    cfg.style = 'both';
    cfg.colormap = custom_cmap;
    cfg.highlight = 'on';
    cfg.highlightchannel = significant_channels{t};
    cfg.highlightsymbol = '.';
    cfg.highlightcolor = 'k';
    cfg.highlightsize = 23;
    cfg.comment = 'no';
    ft_topoplotER(cfg, ds_group_ft);
    pos = get(ax, 'Position');
    pos(3) = 0.14;  
    pos(1) = (t - 1)*0.075; 
    set(ax, 'Position', pos);
    % ft_topoplotER(cfg, ds_z_ft);
end

end
