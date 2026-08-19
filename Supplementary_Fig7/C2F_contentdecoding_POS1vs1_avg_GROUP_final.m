%% MEEG timeseries classification

clear; close all
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/ADAM-1.14-beta/install');
% startup;
%% Load each subject:
numSubjects = 33;
subjectslist = setdiff(1:numSubjects,[22 32 33]); %exclude subject 27
% subjectslist = 1:20;

conditionlist = {'Normal','Masked','LSF'};
figure;
set(gcf,'Color','w');
for cc = 1:length(conditionlist)
clearvars avg group_mvpaoutput sig sigtime
condition = conditionlist{cc};
time_radius = 0;        
trial_bin_num = 4;
nchunks = 10;
OCCIP = 0;
run_pca = 0;

%plot colors
plotcolor= 'b';

counter = 0;
for posvalue1 = 1:4
    for posvalue2 = posvalue1:4
        if posvalue1 ~= posvalue2
        counter = counter + 1;
        for s = subjectslist
            subjNum = s;
            if length(num2str(subjNum)) == 1 || ...
                    (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
                subjNumStr = ['0' num2str(subjNum)];
            else
                subjNumStr = num2str(subjNum);
            end
            
            save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan/PositionDecoding1vs1/',condition,'/noPCA_35Hz_100Hz_', num2str(nchunks), 'chunks_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_position',num2str(posvalue1),num2str(posvalue2),'_LDAba_AllTrials/']);
            if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
            group_mvpaoutput{s,counter} = load([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat'],'sl_map').sl_map;
        
            else
                continue;
            end
        end
    end
    end
end

sub_counter = 0;
for s = 1:size(group_mvpaoutput,1)
    %remove empties
    subjgrp = group_mvpaoutput(s,:);
    subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
    if ~isempty(subjgrp)
        sub_counter = sub_counter + 1;
        submatrix = cell2mat(cellfun(@(x) x.samples, subjgrp, 'UniformOutput', false)');
        avg{1,sub_counter}.samples = mean(submatrix,1);
        avg{1,sub_counter}.a = subjgrp{1}.a;
        avg{1,sub_counter}.fa = subjgrp{1}.fa;
        avg{1,sub_counter}.sa = subjgrp{1}.sa;
        temp_actpattern(:,:,sub_counter) = group_mvpaoutput{s}.activation_pattern;
    end
end

% avg = avg(~cellfun(@isempty, avg));
n_subj = length(avg);

%% Group datasets together:
%some "non-unique" elements due to 0 being jittered, just use the
%parameters from subject 1 because they are mostly the same.

n_subj = length(avg);
ds_cell = cell(n_subj,1);
for p=1:n_subj
        actual_dataset{p,1}  = avg{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        %round values to the 3rd decimal because "0" is a very small number
        %that jitters around 0.
        actual_dataset{p,1}.a.fdim.values{1,1} = round(avg{p}.a.fdim.values{1,1},2);
end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);

ds_group = cosmo_stack(ds_group);

%% Run statistics
%[Optional] - generate a null distribution.
% n_null = 1000;
% null_cell=cell(n_null,1);
% for i_null=1:n_null
%     subj_cell=cell(1,n_subj);
%     for p=1:n_subj
%             ds = group_mvpaoutput{p};
%             ds.sa.targets=1;
%             ds.sa.chunks=p;
%             ds.a.fdim.values{1,1} = round(group_mvpaoutput{p}.a.fdim.values{1,1},3);
%             subj_cell{p}=ds;
%     end
%     null_cell{i_null} = cosmo_stack(subj_cell);
% end
allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group, 'time',allow_clustering_over_time);

% ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',50000);
ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',10000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);

sig = find(ds_z.samples>1.96);

%% Plotting

subplot(3,1,cc)
hold on
timewindow = linspace(-0.4,0.8,241);
% timewindow = -0.4:0.01:0.8;
clearvars temp
plottime = timewindow;
indstart = find(round(avg{1,1}.a.fdim.values{1,1},2) == plottime(1));
indend = find(round(avg{1,1}.a.fdim.values{1,1},2) == plottime(end));

for s = 1:n_subj
    temp(s,:) = avg{s}.samples;
end

meansub = mean(temp);
sem_sub = std(temp)/sqrt(n_subj);
x = timewindow;
curve1 = meansub + sem_sub;
curve2 = meansub - sem_sub;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

h = fill(x2, inBetween,[0 .5 .5]);
set(h, 'EdgeColor', 'none'); % Remove black border
alpha(0.3);

time_values = timewindow;

plot(plottime,meansub,'Color',[0 .5 .5],'LineWidth',1.5);

xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
ylabel('Balanced Accuracy (chance=.5)');

y_pos(1) = 0.465;
y_pos(2) = 0.475;
y_pos(3) = 0.485;

sigtime = time_values(sig);
% onsettime{rr}(condcombo,aa) = time_values(sig(1));
% [~, peakind] = max(meansub);
% peaktime_avg{rr}(condcombo,aa) = time_values(peakind);

if ~isempty(sigtime)
    for i = 1:length(sig)
        % Draw a horizontal line across the full x-axis at each significant value
        line([sigtime(i)-0.005 sigtime(i)+0.005], [y_pos(1) y_pos(1)], ...
             'Color', [0 .5 .5], 'LineWidth', 3);
    end
end
sigmask = abs(ds_z.samples)>1.96;
bold_meansub = meansub;
bold_meansub(~sigmask) = NaN; % set mask to NaN where it is 0
%make significant section bold:
if sum(sigmask) ~= 0; plot(plottime, bold_meansub, 'Color',[0 .5 .5], 'LineWidth', 3.5); end
for i = 1:length(sigtime)
    line([sigtime(i)-0.005 sigtime(i)+0.005], [y_pos(1) y_pos(1)], ...
         'Color', [0 .5 .5], 'LineWidth', 3);
end
% end

xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
ylabel('Balanced accuracy (chance=.5)');
xlabel('time (s)');
title(condition);
xlim([-0.4 0.8]);
% ylim([0.46 0.58]);

if ~isempty(sigtime)
%print first sig. time
x_pos = sigtime(1);
firstsigtime_label = round(sigtime(1)*1000,0); 
text(x_pos, y_pos(1)+0.01, [num2str(firstsigtime_label) 'ms'], 'FontSize', 14, 'Color',[0 .5 .5]);
set(gca,'fontsize', 14)
% ylim([0.48 0.55]);
end

% 
% xlabel('time (s)');
% title(sprintf(['Position decoding', ' - ', condition, ' (N = ', num2str(n_subj),')']));
% xlim([-0.4 0.8]);
% ylim([0.4 0.68]);
% 
% %print first sig. time
% if ~isempty(sig)
%     sigtime = time_values(sig);
%     %make significant section bold:
%     for i = 1:length(sig)
%         % Draw a horizontal line across the full x-axis at each significant value
%         line([sigtime(i)-0.005 sigtime(i)+0.005], [0.48 0.48], ...
%              'Color', [0 .5 .5], 'LineWidth', 3);
%     end
%     x_pos = sigtime(1);
%     y_pos = 0.485;
%     firstsigtime_label = round(sigtime(1)*1000,0); 
%     text(x_pos, y_pos, [num2str(firstsigtime_label) 'ms'], 'FontSize', 10);
% end
% ylim([0.48 0.55]);

% figure; hold on;
% plot(rand(10,1),'b'); 
% plot(rand(10,1),'r');
% plot(rand(10,1),'m');
% plot(rand(10,1),'g');
% legend('Top', 'Bottom', 'Left', 'Right');

% grpmean_actpattern = zscore(grpmean_actpattern, 0, 1); %grpmean_actpattern is 120x241; spatially zscore (across sensors).
% 
% % A = temp_actpattern;
% % temp_actpattern = (A - mean(A, 2)) ./ std(A, [], 2); %z-score
% temp_actpattern = zscore(temp_actpattern, 0, 1); %temp_actpattern is 120x241x30; spatially zscore (across sensors).
% avWeights = squeeze(nanmean(temp_actpattern,3)); % average over subjects
% 
% %% just for the layout labels
% load('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_input/sub33_dataEEG_Rec.mat');
% 
% easycaplayout = load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/layout_acticap128chan.mat').lay_acticap128chan;
% common_labels = intersect(easycaplayout.label,data3.label);
% [~, idx] = ismember(common_labels, easycaplayout.label);
% new_layout = struct();
% 
% new_layout.pos = easycaplayout.pos(idx, :);
% new_layout.width = easycaplayout.width(idx);
% new_layout.height = easycaplayout.height(idx);
% new_layout.label = easycaplayout.label(idx);
% new_layout.outline = easycaplayout.outline;
% new_layout.mask = easycaplayout.mask;
% new_layout.cfg = easycaplayout.cfg;
% layout = new_layout;
% 
% %% Stats for topoplot (activation maps):
% 
% warning('off');
% mpcompcor_method = 'cluster_based';%cluster_based or uncorrected
% indiv_pval = .05;
% cluster_pval = .05;
% 
% clusterPvals = [];
% % load chan_locs in ADAM format:
% % load('/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/channel_info/chanlocs.mat');
% load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/toolboxes/electrodes location files/actiCAP 128 Channel/acticap128.mat');
% eeg_data = load('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_input/sub01_dataEEG_ALL.mat'); %sample data
% label_120elec = eeg_data.data3.label;
% 
% tempchanlocs = acticap128.chanlocs;
% keep_label = ismember({tempchanlocs.labels}, label_120elec);
% chanlocs2 = struct();
% fields = fieldnames(tempchanlocs);  % Get all field names
% for i = 1:numel(fields)
%     temp = {tempchanlocs.(fields{i})};
%     temp2 = temp(keep_label);
%     chanlocs2.(fields{i}) = temp2';
% end
% 
% chanlocs3 = struct('type','EEG','labels',chanlocs2.labels,...
%     'sph_radius',chanlocs2.sph_radius,...
%     'sph_theta',chanlocs2.sph_theta,...
%     'sph_phi',chanlocs2.sph_phi,...
%     'theta',chanlocs2.theta,...
%     'radius',chanlocs2.radius,...
%     'X',chanlocs2.X,'Y',chanlocs2.Y,'Z',chanlocs2.Z);
% 
% chanlocs = chanlocs3;
% settings.chanlocs = chanlocs;
% cfg.tail = 'both'; %'left','both';
% cfg.iterations = 1000;
% 
% % subjweights is a num_subj x channel variable (shold be 30 x 120)
% timeval = linspace(-0.4,0.8,241);
% timelim_all = {[1 80],[81 90],[91 100],[101 110],[111 120],[121 130],[131 140],[141 150],...
%     [151 160],[161 170],[171 180],[181 190],[191 200],[201 210],[211 220],[221 230],[231 240]};
% 
% % timelim_all = {[1 40],[41 80],[81 120],[121 160],[161 200],[201 240]};
% 
% for t = 1:length(timelim_all)
%     subjweights = squeeze(mean(temp_actpattern(:,timelim_all{t}(1):timelim_all{t}(2),:),2))';
%     subjweights(isnan(subjweights)) = 0;
%     if strcmpi(mpcompcor_method, 'cluster_based')
%         % connectivity = get_connected_electrodes({stats.settings.chanlocs(:).labels});
%         connectivity = get_connected_electrodes(chanlocs);
%         [clusterPvals{t}, pStruct{t}] = cluster_based_permutation(subjweights,0,cfg,settings,[],connectivity);
%         % [ clusterPvals, pStruct ] = cluster_based_permutation(data1,data2_or_chance_level,cfg,settings,mask,connectivity)
%     elseif strcmpi(mpcompcor_method, 'uncorrected')
%         [~,clusterPvals{t}] = ttest(subjweights,0,'tail','right');
%     end
% end
% 
% % plot_data = avWeights;
% for t = 1:length(timelim_all)
%     sig_elects{t} = find(clusterPvals{t}<cluster_pval);
% end
% 
% %% plot topoplots
% 
% %smooth here for visualization:
% % avWeights = imgaussfilt(avWeights, 1);
% 
% timeval = linspace(-0.4,0.8,241);
% timelim_all = {[1 80],[81 90],[91 100],[101 110],[111 120],[121 130],[131 140],[141 150],...
%     [151 160],[161 170],[171 180],[181 190],[191 200],[201 210],[211 220],[221 230],[231 240]};
% % timelim_all = {[1 40],[41 80],[81 120],[121 160],[161 200],[201 240]};
% 
% absval = abs(avWeights);
% maxabsval = max(max(avWeights));
% 
% figure; set(gcf,'Color','w');
% for t = 1:length(timelim_all)
%     activation_temp = nanmean(avWeights(:,timelim_all{t}(1):timelim_all{t}(2)),2);
% 
%     data = [];
%     data.avg = activation_temp;
%     data.label = data3.label;
%     data.dimord = 'chan_time';
%     data.time = 0; 
% 
%     cfg = [];
%     cfg.layout = layout;
%     cfg.marker = 'off';
%     cfg.comment = 'no';
%     cfg.zlim = [-maxabsval maxabsval]; %
%     % cfg.style = 'straight';
%     cfg.highlight = 'on';
%     cfg.highlightchannel = sig_elects{t};  % List of sig. channels to highlight
%     cfg.highlightsymbol = '.';  % Circle marker
%     cfg.highlightstyle = 'marker';
%     cfg.highlightfontsize = 24;
%     subplot(4,5,t);
%     ft_topoplotER(cfg, data);
%     h = cbar('vert');
%     bartitle = 'spatially z-scored';
%     ylabel(h,bartitle);
%     cmap  = brewermap([],'*RdBu');
%     colormap(gcf,cmap);
%     title(sprintf([num2str(round(timeval(timelim_all{t}(1)),2)),' to ', num2str(round(timeval(timelim_all{t}(2)),2)), 's']));
% end
% nosedir = '+X';
% for t = 1:length(timelim_all)
%     subplot(4,5,t);
%     elecs = sig_elects{t};
%     activation_temp = mean(grpmean_actpattern(:,timelim_all{t}(1):timelim_all{t}(2)),2);
%     weightlim = [-maxabsval maxabsval];
%     topoplot_jjf(activation_temp,convertlocs(chanlocs,'cart2topo'),'maplimits',weightlim,'style','map','electrodes','on','plotrad',.65,'nosedir',nosedir,'emarker2',{elecs,'o','k',5,1}); % 'electrodes','ptslabels', 'plotrad',.7
%     title(sprintf([num2str(round(timeval(timelim_all{t}(1)),2)),' to ', num2str(round(timeval(timelim_all{t}(2)),2)), 's']));
%     h = cbar('vert');
%     bartitle = 'spatially z-scored';
%     ylabel(h,bartitle);
%     cmap  = brewermap([],'*RdBu');
%     colormap(gcf,cmap);
% end

% %% plot topoplots
% 
% %smooth here for visualization:
% grpmean_actpattern = imgaussfilt(grpmean_actpattern, 1);
% 
% timeval = linspace(-0.4,0.8,241);
% timelim_all = {[1 40],[51 60],[61 70],[71 80],[81 90],[91 100],[101 110],[111 120],[121 130],[131 140],[141 150],...
%     [151 160],[161 170],[171 180],[181 190],[191 200],[201 210],[211 220],[221 230],[231 240]};
% 
% absval = abs(grpmean_actpattern);
% maxabsval = max(max(absval));
% 
% figure; set(gcf,'Color','w'); hold on
% for t = 1:length(timelim_all)
% activation_temp = mean(grpmean_actpattern(:,timelim_all{t}(1):timelim_all{t}(2)),2);
% 
% data = [];
% data.avg = activation_temp;
% data.label = data3.label;
% data.dimord = 'chan_time';
% data.time = 0;  % Single timepoint
% 
% % Step 3: Plot using ft_topoplotER
% cfg = [];
% cfg.layout = layout;
% cfg.marker = 'on';
% cfg.comment = 'no';
% cfg.zlim = [-maxabsval maxabsval]; %
% cfg.style = 'straight';
% subplot(4,5,t);
% ft_topoplotER(cfg, data);
% title(sprintf([num2str(round(timeval(timelim_all{t}(1)),2)),' to ', num2str(round(timeval(timelim_all{t}(2)),2)), 's']));
% 
% cmap  = brewermap([],'*RdBu');
% colorbar();
% colormap(gcf,cmap);
% end
end
