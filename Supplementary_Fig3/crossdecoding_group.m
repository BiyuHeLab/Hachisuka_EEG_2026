%% Cross-condition decoding: Group Analysis

% Train on Normal/Masked/LSF, and test on Normal/Masked/LSF
% For both superordinate and basic-level decoding.
% Plots group-level analysis and bootstrapped onset & peak times

% Subject-level analysis: C2F_crossdecoding_leaveoutexemp_batch.m
% Bootstrapping: C2F_bootstrap_peakonset_groupfinal_crosscond.m

% Last updated July 1 2025, Ayaka Hachisuka (ahachisu@gmail.com)

clear;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');

numSubjects = 33;
subjectslist = setdiff(1:numSubjects,[22 32 33]);

run_pca = 0;
reverse = 1;
superorbasicLIST = {'Basic','Super'};
trialtype = 'ALL'; %'rec' or 'unrec'

n_modalities = 2; %it's a 2x2 cross-decoding schematic.

exemppairlist{1} = 1:3;
exemppairlist{2} = 4:6;
exemppairlist{3} = 7:9;
exemppairlist{4} = 10:12;

conditionlist1{1} = 'Normal'; conditionlist2{1} = 'LSF';
conditionlist1{2} = 'Masked'; conditionlist2{2} = 'LSF';
conditionlist1{3} = 'Normal'; conditionlist2{3} = 'Masked';

figure; set(gcf,'Color','w');
hold on

figcounter = 0;

for rr = 1:2
    reverse = rr - 1;
    % plot1 = gobjects(1,3);
for condcombo = 1:length(conditionlist1)
    figcounter = figcounter + 1;
    condition1 = conditionlist1{condcombo};
    condition2 = conditionlist2{condcombo};
for aa = 1:length(superorbasicLIST)
    superorbasic = superorbasicLIST{aa};
    counter1 = 0; counter2 = 0;
for s = subjectslist
    for pair = 1
        for train_modality=1:n_modalities
            for test_modality=1:n_modalities
                if train_modality ~= test_modality
                   if train_modality > test_modality && reverse == 0
                       % if train_modality > test_modality %if reverse
                        counter2 = counter2 + 1;
                        if contains(superorbasic,'Exemplar')
                            exempcounter = 0;
                            for pp = 1:4
                            exemp_pair = exemppairlist{pp};
                            for exempval1 = exemp_pair(1):exemp_pair(3) %one vs other schematic
                                for exempval2 = exempval1:exemp_pair(3) %one vs other schematic
                                    if exempval1 ~= exempval2
                                        exempcounter = exempcounter + 1;
                            sl_result_group_exemp{exempcounter}{counter2} = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'/',condition1,'_',condition2,'sub',num2str(s),'_', ...
                num2str(train_modality), num2str(test_modality), '_exemp',num2str(exempval1),num2str(exempval2),'_ds_searchlight_result']),'ds_searchlight_result');
                                    end
                                end
                            end
                            end
                        elseif contains(superorbasic,'Basic')
                            sl_result_group1a{counter2} = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'1/',condition1,'_',condition2,'sub',num2str(s),'_', ...
                num2str(train_modality), num2str(test_modality), '_300repeats_ds_searchlight_result']),'ds_searchlight_result');
                            sl_result_group1b{counter2} = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'2/',condition1,'_',condition2,'sub',num2str(s), '_',...
                num2str(train_modality), num2str(test_modality), '_300repeats_ds_searchlight_result']),'ds_searchlight_result');
                        elseif contains(superorbasic,'Super')
                        sl_result_group1{counter2} = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'/',condition1,'_',condition2,'sub',num2str(s), '_',...
                num2str(train_modality), num2str(test_modality), '_300repeats_ds_searchlight_result']),'ds_searchlight_result');
                        end
                   elseif train_modality > test_modality && reverse == 1
                       % if train_modality > test_modality %if reverse
                        counter2 = counter2 + 1;
                        if contains(superorbasic,'Exemplar')
                            exempcounter = 0;
                            for pp = 1:4
                            exemp_pair = exemppairlist{pp};
                            exemp_counter = 0;
                            for exempval1 = exemp_pair(1):exemp_pair(3) %one vs other schematic
                                for exempval2 = exempval1:exemp_pair(3) %one vs other schematic
                                    if exempval1 ~= exempval2
                                        exempcounter = exempcounter + 1;
                            sl_result_group_exemp{exempcounter}{counter2} = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'/',condition1,'_',condition2,'sub',num2str(s),'_', ...
                num2str(train_modality), num2str(test_modality), '_exemp',num2str(exempval1),num2str(exempval2),'_ds_searchlight_result']),'ds_searchlight_result');
                                    end
                                end
                            end
                            end
                        elseif contains(superorbasic,'Basic')
                            sl_result_group1a{counter2} = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'1/',condition1,'_',condition2,'sub',num2str(s),'_', ...
                num2str(test_modality), num2str(train_modality), '_300repeats_ds_searchlight_result']),'ds_searchlight_result');
                            sl_result_group1b{counter2} = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'2/',condition1,'_',condition2,'sub',num2str(s), '_',...
                num2str(test_modality), num2str(train_modality), '_300repeats_ds_searchlight_result']),'ds_searchlight_result');
                        elseif contains(superorbasic,'Super')
                        sl_result_group1{counter2} = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/FinalCrossDecode_results/LDA/AllChan200Hz_balexemp/Cross_Cond_decoding_LDAba_',trialtype,'/',...
                superorbasic,'/',condition1,'_',condition2,'sub',num2str(s), '_',...
                num2str(test_modality), num2str(train_modality), '_300repeats_ds_searchlight_result']),'ds_searchlight_result');
                        end
                    end
                end
            end
        end
    end
end

%% Average within-decoding subjects:
%if it's for BASIC category decoding, combine the two 
counter1 = 0;
counter2 = 0;
counter3 = 0;
if contains(superorbasic,'Basic')
    for s = 1:size(sl_result_group1a,2)
        %remove empties
        subjgrp = sl_result_group1a(s);
        subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
        if ~isempty(subjgrp)
            counter1 = counter1 + 1;
            avg1a{1,counter1}.samples = subjgrp{1,1}.ds_searchlight_result.samples;
            % avg1a{1,counter1}.a = subjgrp{1,1}.ds_searchlight_result.a;
            % avg1a{1,counter1}.fa = subjgrp{1,1}.ds_searchlight_result.fa;
            % avg1a{1,counter1}.sa = subjgrp{1,1}.ds_searchlight_result.sa;
        end
    end
    for s = 1:size(sl_result_group1b,2)
        subjgrp = sl_result_group1b(s);
        subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
        if ~isempty(subjgrp)
            counter2 = counter2 + 1;
            avg1b{1,counter2}.samples = subjgrp{1,1}.ds_searchlight_result.samples;
            % avg1a{1,counter2}.a = subjgrp{1,1}.ds_searchlight_result.a;
            % avg1a{1,counter2}.fa = subjgrp{1,1}.ds_searchlight_result.fa;
            % avg1a{1,counter2}.sa = subjgrp{1,1}.ds_searchlight_result.sa;
        end
    end
    for s = 1:size(sl_result_group1a,2)
        newsubjgrp{1,s}.ds_searchlight_result.samples = nanmean([avg1a{1,s}.samples; avg1b{1,s}.samples]);
        avg1{1,s} = newsubjgrp{1,s}.ds_searchlight_result;
    end
elseif contains(superorbasic,'Exemplar')
    for s = 1:size(sl_result_group_exemp{1},2)
        for exempcounter = 1:size(sl_result_group_exemp,2)
            %remove empties
            subjgrp = sl_result_group_exemp{exempcounter}(s);
            subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
            if ~isempty(subjgrp)
                allexempsamples(exempcounter,:) = subjgrp{1,1}.ds_searchlight_result.samples;
            end
        end
        avg1a{1,s}.samples = nanmean(allexempsamples);
    end
   
    for s = 1:size(sl_result_group_exemp{1},2)
        newsubjgrp{1,s}.ds_searchlight_result.samples = avg1a{1,s}.samples;
        avg1{1,s} = newsubjgrp{1,s}.ds_searchlight_result;
    end
else
    counter = 0;
    for s = 1:size(sl_result_group1,2)
        %remove empties
        subjgrp = sl_result_group1(s);
        subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
        if ~isempty(subjgrp)
            counter = counter + 1;
            avg1{1,counter}.samples = subjgrp{1,1}.ds_searchlight_result.samples;
            % avg1{1,counter}.a = subjgrp{1,1}.ds_searchlight_result.a;
            % avg1{1,counter}.fa = subjgrp{1,1}.ds_searchlight_result.fa;
            % avg1{1,counter}.sa = subjgrp{1,1}.ds_searchlight_result.sa;
        end
    end
end

%% Group datasets together:
%some "non-unique" elements due to 0 being jittered, just use the
%parameters from subject 1 because they are mostly the same.

%within
clearvars actual_dataset
counter = 0;

for p=1:length(avg1)
    actual_dataset{p,1}.samples = avg1{p}.samples;
    actual_dataset{p,1}.sa.targets = 1;
    actual_dataset{p,1}.sa.chunks = p;
    actual_dataset{p,1}.fa.time = (1:241);
    actual_dataset{p,1}.fa.center_ids = (1:241);
    actual_dataset{p,1}.a.meeg.samples_field = 'trial';
    actual_dataset{p,1}.a.fdim.labels = {'time'};
    actual_dataset{p,1}.a.fdim.values{1,1} = linspace(-0.4,0.8,241);
end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group_within = cosmo_stack(ds_group);

allbalacc{condcombo,aa} = ds_group_within.samples;

allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group_within, 'time',allow_clustering_over_time);
% ds_z = cosmo_montecarlo_cluster_stat(ds_group_within,nbrhood,'h0_mean',0.5,'niter',1000);
ds_z = cosmo_montecarlo_cluster_stat(ds_group_within,nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);
sig_within = find(abs(ds_z.samples)>1.96);
sigmask = abs(ds_z.samples)>1.96;
%% PLOTTING
subplot(2,3,figcounter); hold on;
timewindow = linspace(-0.4,0.8,241); %-0.4:0.01:0.8;
clearvars temp
plottime = timewindow;
indstart = find(round(actual_dataset{1,1}.a.fdim.values{1,1},3) == plottime(1));
indend = find(round(actual_dataset{1,1}.a.fdim.values{1,1},3) == plottime(end));

for s = 1:length(subjectslist)
    temp(s,:) = avg1{s}.samples;%(indstart(1):indend(end));
end

meansub = nanmean(temp);
sem_sub = nanstd(temp)/sqrt(length(subjectslist));
x = timewindow;
curve1 = meansub + sem_sub;
curve2 = meansub - sem_sub;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];

if contains(superorbasic,'Super')
    plotcolor = [50 130 246]./255;
elseif contains(superorbasic,'Basic')
    plotcolor = [240 134 80]./255;
elseif contains(superorbasic,'Exemplar')
    plotcolor = [58 224 213]./255;
end


h = fill(x2, inBetween, plotcolor);
set(h, 'EdgeColor', 'none'); % Remove black border
alpha(0.15);

time_values = timewindow;
plot1(aa) = plot(plottime,meansub,'Color',plotcolor,'LineWidth',1);

%make significant section bold:
if sum(sigmask) ~= 0
    boldvals = meansub; 
    boldvals(~sigmask) = NaN; 
    plot(plottime, boldvals, 'Color',plotcolor, 'LineWidth', 3.5); 
end

% y_pos(1) = 0.465;
y_pos(1) = 0.485;
y_pos(2) = 0.495;

sigtime = time_values(sig_within);
if ~isempty(sigtime)
    for c = 1:length(sig_within)
        % Draw a horizontal line across the full x-axis at each significant value
        line([sigtime(c)-0.0025 sigtime(c)+0.0025], [y_pos(aa) y_pos(aa)], ...
             'Color', plotcolor, 'LineWidth', 3);
    end
end

xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
ylabel('Balanced accuracy (chance=.5)');
xlabel('time (s)');
grid 'on'
xlim([-0.4 0.8]);
ylim([0.48 0.56]);

if ~isempty(sigtime)
%print first sig. time
x_pos = sigtime(1);
firstsigtime_label = round(sigtime(1)*1000,0); 
text(x_pos, y_pos(aa)+0.0025, [num2str(firstsigtime_label) 'ms'], 'FontSize', 10, 'Color',plotcolor);
% ylim([0.48 0.55]);
end

end

diff_acc = allbalacc{condcombo,1} - allbalacc{condcombo,2};
ds_group_within.samples = diff_acc;
allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group_within, 'time',allow_clustering_over_time);
% ds_z = cosmo_montecarlo_cluster_stat(ds_group_within,nbrhood,'h0_mean',0,'niter',1000);
ds_z = cosmo_montecarlo_cluster_stat(ds_group_within,nbrhood,'h0_mean',0.5,'niter', 1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);

sig_ind = find(abs(ds_z.samples)>1.96);
% 
% thresh = 0.05;
% for t = 1:241
%     [binary(t),pval(t)] = ttest(allbalacc{condcombo,1}(:,t),allbalacc{condcombo,2}(:,t));
% end
% orig_binary = double(pval < thresh);
% [clust_size,start,stop,clust_labels] = find_clusters(orig_binary);
% num_perm = 10000;
% for p = 1:num_perm
%     shuffled_binary = orig_binary(randperm(length(orig_binary)));
%     [clust_size_shuffled,~,~,~] = find_clusters(shuffled_binary);
%     permuted_max_cluster(p) = max(clust_size_shuffled);
% end
% 
% sig_cluster = size(clust_size,1);
% for clustind = 1:length(clust_size)
%     clust_stats_pval = sum(permuted_max_cluster > clust_size(clustind)) / num_perm;
%     if clust_stats_pval < 0.025
%         sig_cluster(clustind) = 1; 
%         sig_clust_start(clustind) = start(clustind);
%         sig_clust_end(clustind) = stop(clustind);
%     else
%         sig_cluster(clustind) = 0; 
%         sig_clust_start(clustind) = NaN;
%         sig_clust_end(clustind) = NaN;
%     end
% end
% 
% sig_ind = arrayfun(@(s, e) s:e, sig_clust_start, sig_clust_end, 'UniformOutput', false);
% sig_ind = [sig_ind{:}];  % Concatenate into a single vector
% sig_ind = sig_ind(~isnan(sig_ind));
% timeval = linspace(-0.4,0.8,241);

corrsigtimeall = plottime(sig_ind);
if ~isempty(corrsigtimeall)
for c = 1:length(corrsigtimeall)
line([corrsigtimeall(c)-0.0025 corrsigtimeall(c)+0.0025], [0.54 0.54], ...
         'Color', 'k', 'LineWidth', 3);
end

diffsig_label = round(corrsigtimeall(1)*1000,0); 
x_pos = plottime(sig_ind(1));
text(x_pos, 0.544, [num2str(diffsig_label) 'ms'], 'FontSize', 10);

end

num_subjects = length(sl_result_group1);
if reverse == 0
    % legend([plot1(1),plot1(2)],{sprintf([condition1,'-',condition2, ' (Super)']), sprintf([condition1,'-',condition2, ' (Basic)'])});
    title(sprintf([condition1, ' - ', condition2, ' (N = ',num2str(num_subjects),')']));
else
    % legend([plot1(1),plot1(2)],{sprintf([condition2,'-',condition1, ' (Super)']), sprintf([condition2,'-',condition1, ' (Basic)'])});
    title(sprintf([condition2, ' - ', condition1, ' (N = ',num2str(num_subjects),')']));
end
end
end

%% Bootstrapped onsets:
bootstrapNUM = 1000;
additionalCOND = [];
condColor{1} = [50 130 246]./255;
condColor{2} = [240 134 80]./255;
condColor{3} = [58 224 213]./255;
conditionList1 = {'Normal','Normal','Masked'};
conditionList2 = {'Masked','LSF','LSF'};
superorbasicLIST = {'Super','Basic','Exemplar'};
trialtype = 'ALL';
time_radius = 0;        
trial_bin_num = 4;
numrepeats = 300;
OCCIP = 0;
run_pca = 0;
n_modalities = 2;
timewindow = linspace(-0.4,0.8,241);

clearvars confidence_int
for reverse = 1:2
figure; set(gcf,'Color','w'); hold on

for aa = 1:length(superorbasicLIST)
analysistype = superorbasicLIST{aa};
for c = 1:3
    
    cond1 = conditionList1{c};
    cond2 =  conditionList2{c};
    
    if reverse ==1 
        savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/CrossCond_Bootstrap/Bootstrapping1000_LDA_',trialtype,'Trials_stats/',cond1,'_',cond2,'/'];
    elseif reverse == 2
        savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/CrossCond_Bootstrap/Bootstrapping1000_LDA_',trialtype,'Trials_stats/',cond2,'_',cond1,'/'];
    end


    temp = load(fullfile([savedir,'ds_z_',analysistype]),'ds_z');
    ds_z_all = temp.ds_z;
    for ii = 1:length(ds_z_all)
        sigvals = find(abs(ds_z_all{ii}.samples)>1.96);
        if ~isempty(sigvals)
            sig_onset(ii) = timewindow(sigvals(1));
        else
            sig_onset(ii) = nan;
        end
    end

    onsettime(c,aa) = nanmean(sig_onset);
    sig_onset_all{c,aa} = sig_onset;

    confidence_int(1,aa) = prctile(sig_onset, 2.5); %low
    confidence_int(2,aa) = prctile(sig_onset, 97.5); %high
  
    CIbyCond_onset{c} = confidence_int; %CItdist; %confidence_int %yCI95
end

end


% yMeanbyCond_onset{c} = yMean;

for c = 1:3
subplot(1,3,c); hold on
bar(1,onsettime(c,1),'FaceColor',condColor{1});
bar(2,onsettime(c,2),'FaceColor',condColor{2})
bar(3,onsettime(c,3),'FaceColor',condColor{3})
errorbar(1,onsettime(c,1), onsettime(c,1) - CIbyCond_onset{c}(1,1),CIbyCond_onset{c}(2,1) - onsettime(c,1),'.','Color','k','Linewidth',1);
errorbar(2,onsettime(c,2), onsettime(c,2) - CIbyCond_onset{c}(1,2),CIbyCond_onset{c}(2,2) - onsettime(c,2),'.','Color','k','Linewidth',1);
errorbar(3,onsettime(c,3), onsettime(c,3) - CIbyCond_onset{c}(1,3),CIbyCond_onset{c}(2,3) - onsettime(c,3),'.','Color','k','Linewidth',1);
xlim([0 5]);
ylim([-0.4 0.8])
xticks(1:3);
xticklabels({'Super','Basic','Exemplar'});
xtickangle(90);
ylabel('Onset Latency (ms)');
if reverse == 1
    % legend([plot1(1),plot1(2)],{sprintf([condition1,'-',condition2, ' (Super)']), sprintf([condition1,'-',condition2, ' (Basic)'])});
    title(sprintf([cond1, ' - ', cond2, ' (N = ',num2str(numSubjects),')']));
else
    % legend([plot1(1),plot1(2)],{sprintf([condition2,'-',condition1, ' (Super)']), sprintf([condition2,'-',condition1, ' (Basic)'])});
    title(sprintf([cond2, ' - ', cond1, ' (N = ',num2str(numSubjects),')']));
end

for aa = 1:length(superorbasicLIST)
    fprintf([superorbasicLIST{aa}, ': ', cond1, ' - ', cond2,'\n']);
    fprintf('%.4f [%.4f %.4f]\n', onsettime(c,aa), CIbyCond_onset{c}(1,aa), CIbyCond_onset{c}(2,aa));
end
end
end



%% Peak latency test:

for reverse = 1:2
figure; set(gcf,'Color','w'); hold on

for aa = 1:length(superorbasicLIST)
analysistype = superorbasicLIST{aa};
for c = 1:3
    
    cond1 = conditionList1{c};
    cond2 =  conditionList2{c};
    
    if reverse ==1 
        savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/CrossCond_Bootstrap/Bootstrapping1000_LDA_',trialtype,'Trials_stats/',cond1,'_',cond2,'/'];
    elseif reverse == 2
        savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/CrossCond_Bootstrap/Bootstrapping1000_LDA_',trialtype,'Trials_stats/',cond2,'_',cond1,'/'];
    end

    temp = load(fullfile([savedir,'bootstrap_balacc_data_',analysistype]),'bootstrap_data_all');
    ds_z_all = temp.bootstrap_data_all;

    for b = 1:bootstrapNUM
        [~,peakval(b)] = max(nanmean(ds_z_all{b}.samples));
    end
    peaktimes = timewindow(peakval);
    peak_all(c,aa) = nanmean(peaktimes);

    confidence_int(1,aa) = prctile(peaktimes, 2.5); %low
    confidence_int(2,aa) = prctile(peaktimes, 97.5); %high
  
    CIbyCond_peak{c} = confidence_int; %CItdist; %confidence_int %yCI95
end

end


% yMeanbyCond_peak{c} = yMean;

for c = 1:3
subplot(1,3,c); hold on
bar(1,peak_all(c,1),'FaceColor',condColor{1});
bar(2,peak_all(c,2),'FaceColor',condColor{2})
bar(3,peak_all(c,3),'FaceColor',condColor{3})
errorbar(1,peak_all(c,1), peak_all(c,1) - CIbyCond_peak{c}(1,1),CIbyCond_peak{c}(2,1) - peak_all(c,1),'.','Color','k','Linewidth',1);
errorbar(2,peak_all(c,2), peak_all(c,2) - CIbyCond_peak{c}(1,2),CIbyCond_peak{c}(2,2) - peak_all(c,2),'.','Color','k','Linewidth',1);
errorbar(3,peak_all(c,3), peak_all(c,3) - CIbyCond_peak{c}(1,3),CIbyCond_peak{c}(2,3) - peak_all(c,3),'.','Color','k','Linewidth',1);
xlim([0 5]);
ylim([-0.4 0.8])
xticks(1:3);
xticklabels({'Super','Basic','Exemplar'});
xtickangle(90);
ylabel('peak Latency (ms)');
if reverse == 1
    % legend([plot1(1),plot1(2)],{sprintf([condition1,'-',condition2, ' (Super)']), sprintf([condition1,'-',condition2, ' (Basic)'])});
    title(sprintf([cond1, ' - ', cond2, ' (N = ',num2str(numSubjects),')']));
else
    % legend([plot1(1),plot1(2)],{sprintf([condition2,'-',condition1, ' (Super)']), sprintf([condition2,'-',condition1, ' (Basic)'])});
    title(sprintf([cond2, ' - ', cond1, ' (N = ',num2str(numSubjects),')']));
end

for aa = 1:length(superorbasicLIST)
    fprintf([superorbasicLIST{aa}, ': ', cond1, ' - ', cond2,'\n']);
    fprintf('%.4f [%.4f %.4f]\n', peak_all(c,aa), CIbyCond_peak{c}(1,aa), CIbyCond_peak{c}(2,aa));
end
end
end

