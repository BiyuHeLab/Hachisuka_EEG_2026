%% Category decoding cross position group analysis

clear; close all
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/ADAM-1.14-beta/install');
% startup;
%% Load each subject:
numSubjects = 33;
subjectslist = setdiff(1:numSubjects,[22 32 33]); %exclude subject 27
% subjectslist = 1:20;

conditionlist = {'Normal','Masked','LSF'}; %{'AllCond'}; %
% conditionlist = {'AllCond'};
superorbasicLIST = {'Super','Basic'};
% analysistype = 'CategoryDecodingByPos_leaveonetrialout';
analysistype = 'CategoryDecodingCrossPos_leaveonetrialout_300repeats';
figure; set(gcf,'Color','w');
hold on
for cattype = 1:length(superorbasicLIST)

superorbasic = superorbasicLIST{cattype};

cc_counter = 0;

for cc = 1:length(conditionlist)
clearvars avg_vals group_mvpaoutput sig sigtime
condition = conditionlist{cc};
time_radius = 0;        
trial_bin_num = 4;
nchunks = 10;
OCCIP = 0;
run_pca = 0;

%plot colors
plotcolor= 'b';


for posvalue1 = 1:4
        scounter = 0;
        for s = subjectslist
            subjNum = s;
            if length(num2str(subjNum)) == 1 || ...
                    (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
                subjNumStr = ['0' num2str(subjNum)];
            else
                subjNumStr = num2str(subjNum);
            end
           
            if contains(superorbasic,'Basic')
                save_dir1 = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Final_Resubmission_results_April16/Fig3/',analysistype,'/',condition,'/Basic1/35Hz_0timerad_',num2str(trial_bin_num),'_position',num2str(posvalue1),'_LDAba_AllTrials/']);
                save_dir2 = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Final_Resubmission_results_April16/Fig3/',analysistype,'/',condition,'/Basic2/35Hz_0timerad_',num2str(trial_bin_num),'_position',num2str(posvalue1),'_LDAba_AllTrials/']);
                % save_dir1=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan/',analysistype,'/',condition,'/Basic1/noPCA_35Hz_100Hz_', num2str(nchunks), 'chunks_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_position',num2str(posvalue1),'_LDAba_AllTrials/']);
                % save_dir2=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan/',analysistype,'/',condition,'/Basic2/noPCA_35Hz_100Hz_', num2str(nchunks), 'chunks_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_position',num2str(posvalue1),'_LDAba_AllTrials/']);
                if exist(fullfile([save_dir1,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
                    temp1 = load([save_dir1,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat'],'sl_map').sl_map;
                else temp1 = [];
                end
                if exist(fullfile([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
                    temp2 = load([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat'],'sl_map').sl_map;
                else temp2 = [];
                end
                if isempty(temp1) && ~isempty(temp2)
                    scounter = scounter + 1;
                    group_mvpaoutput{scounter,posvalue1} = [temp2.samples];
                elseif isempty(temp2) && ~isempty(temp1)
                    scounter = scounter + 1;
                    group_mvpaoutput{scounter,posvalue1} = [temp1.samples];
                elseif ~isempty(temp2) && ~isempty(temp1)
                    scounter = scounter + 1;
                    group_mvpaoutput{scounter,posvalue1} = nanmean([temp1.samples;temp2.samples]);
                end
            else
                % save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan/',analysistype,'/',condition,'/',superorbasic,'/noPCA_35Hz_100Hz_', num2str(nchunks), 'chunks_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_position',num2str(posvalue1),'_LDAba_AllTrials/']);
                save_dir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Final_Resubmission_results_April16/Fig3/',analysistype,'/',condition,'/Super/35Hz_0timerad_',num2str(trial_bin_num),'_position',num2str(posvalue1),'_LDAba_AllTrials/']);
                if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
                    scounter = scounter + 1;
                    group_mvpaoutput{scounter,posvalue1} = load([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat'],'sl_map').sl_map;
                end
            end
        end
end

clearvars avg_vals
sub_counter = 0;
for s = 1:size(group_mvpaoutput,1)
    %remove empties
    subjgrp = group_mvpaoutput(s,:);
    subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
    if ~isempty(subjgrp)
        sub_counter = sub_counter + 1;
        for p = 1:4
            if contains(superorbasic,'Basic')
                submatrix = cell2mat(cellfun(@(x) x, subjgrp, 'UniformOutput', false)');
            else
                submatrix = cell2mat(cellfun(@(x) x.samples, subjgrp, 'UniformOutput', false)');
            end
        end

        if size(submatrix,1) > 1
            avg_vals{sub_counter}.samples = nanmean(submatrix);
        else
            avg_vals{sub_counter}.samples = submatrix;
        end
        avg_vals{sub_counter}.a.fdim.labels = {'time'};
        avg_vals{sub_counter}.a.fdim.values{1,1} = linspace(-0.4,0.8,241);
        avg_vals{sub_counter}.fa.time = 1:241;
        avg_vals{sub_counter}.fa.center_ids = 1:241;
        avg_vals{sub_counter}.sa.labels = {'balanced_accuracy'};
        % temp_actpattern(:,:,sub_counter) = group_mvpaoutput{s}.activation_pattern;
        
    end
end

% avg = avg(~cellfun(@isempty, avg));
n_subj = length(avg_vals);

%% Group datasets together:
%some "non-unique" elements due to 0 being jittered, just use the
%parameters from subject 1 because they are mostly the same.

n_subj = length(avg_vals);
ds_cell = cell(n_subj,1);

    for p=1:n_subj
        actual_dataset{p,1} = avg_vals{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        %round values to the 3rd decimal because "0" is a very small number
        %that jitters around 0.
        actual_dataset{p,1}.a.fdim.values{1,1} = linspace(-0.4,0.8,241);
    end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);

ds_group = cosmo_stack(ds_group);

allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group, 'time',allow_clustering_over_time);

% ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',50000);
ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',10000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);

sig = find(ds_z.samples>1.96);

%% Plotting

if contains(superorbasic,'Basic')
    plotcolor = [240 134 80]./255;
    y_pos(1) = 0.465;
else
    plotcolor = [50 130 246]./255;
    y_pos(1) = 0.485;
end

cc_counter = cc_counter + 1;
if ~strcmp(condition,'AllCond')
    subplot(3,1,cc_counter)
end
hold on
timewindow = linspace(-0.4,0.8,241);
% timewindow = -0.4:0.01:0.8;
clearvars temp
plottime = timewindow;
indstart = find(round(avg_vals{1}.a.fdim.values{1,1},2) == plottime(1));
indend = find(round(avg_vals{1}.a.fdim.values{1,1},2) == plottime(end));

for s = 1:n_subj
    temp(s,:) = avg_vals{s}.samples;
end

meansub = nanmean(temp);
sem_sub = nanstd(temp)/sqrt(n_subj);
x = timewindow;
curve1 = meansub + sem_sub;
curve2 = meansub - sem_sub;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

h = fill(x2, inBetween,plotcolor);
set(h, 'EdgeColor', 'none'); % Remove black border
alpha(0.3);

time_values = timewindow;

plot(plottime,meansub,'Color',plotcolor,'LineWidth',1);

% if pos == 1
%     position = 'Top';
% elseif pos == 2
%     position = 'Bottom';
% elseif pos == 3
%     position = 'Left';
% elseif pos == 4
%     position = 'Right';
% end
% title(sprintf([condition, '-', position]));
title(sprintf([condition]));
xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
ylabel('Balanced Accuracy (chance=.5)');

sigtime = time_values(sig);
% onsettime{rr}(condcombo,aa) = time_values(sig(1));
% [~, peakind] = max(meansub);
% peaktime_avg_vals{rr}(condcombo,aa) = time_values(peakind);

if ~isempty(sigtime)
    for i = 1:length(sig)
        % Draw a horizontal line across the full x-axis at each significant value
        line([sigtime(i)-0.005 sigtime(i)+0.005], [y_pos(1) y_pos(1)], ...
             'Color', plotcolor, 'LineWidth', 3);
    end
end
sigmask = abs(ds_z.samples)>1.96;
bold_meansub = meansub;
bold_meansub(~sigmask) = NaN; % set mask to NaN where it is 0
%make significant section bold:
if sum(sigmask) ~= 0; plot(plottime, bold_meansub, 'Color',plotcolor, 'LineWidth', 3.5); end
for i = 1:length(sigtime)
    line([sigtime(i)-0.005 sigtime(i)+0.005], [y_pos(1) y_pos(1)], ...
         'Color', plotcolor, 'LineWidth', 3);
end
% end

xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
ylabel('Balanced accuracy (chance=.5)');
xlabel('time (s)');

xlim([-0.4 0.8]);
ylim([0.46 0.545]);
grid 'on';
if ~isempty(sigtime)
%print first sig. time
x_pos = sigtime(1);
firstsigtime_label = round(sigtime(1)*1000,0); 
text(x_pos, y_pos(1)+0.005, [num2str(firstsigtime_label) 'ms'], 'FontSize', 10, 'Color',plotcolor);
end

ax = gca;
ax.XTick = -0.4:0.2:0.8;
ax.XMinorTick = 'on';
ax.XAxis.MinorTickValues = -0.4:0.1:0.8;

set(gca,'Fontsize',14);
end
end
