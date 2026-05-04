%% Calculating group-level analyses for category decoding:

% Last updated Feb 6 2026, Ayaka Hachisuka

clear; close all
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
%% Load each subject:
numSubjects = 33;
conditionList = {'Normal','Masked','LSF'};
timewindow = linspace(-0.4,0.8,241);
plotcounter = 0;

figure; set(gcf,'Color','w');
for c = 1:length(conditionList)
condition =  conditionList{c};

% For plotting bar charts for onset & peak times: include exemplar
% decoding:
% superorbasicLIST = {'Super','Basic','Exemplar'};
% analysistypeLIST = {'Super','Basic','Exemplar'};
superorbasicLIST = {'Super','Basic'};
analysistypeLIST = {'Super','Basic'};

trialtype = 'ALL';
time_radius = 0;        
trial_bin_num = 4;
numrepeats = 300;
OCCIP = 0;
run_pca = 0;
subjectslist = setdiff(1:numSubjects,[22 32 33]); 

exemppairlist{1} = 1:3;
exemppairlist{2} = 4:6;
exemppairlist{3} = 7:9;
exemppairlist{4} = 10:12;

subplot(3,1,c);
hold on
for aa = 1:length(superorbasicLIST)
    clearvars ds_group
    superorbasic = superorbasicLIST{aa};
    analysistype = analysistypeLIST{aa};
    clearvars group_mvpaoutput
    for s = subjectslist
    
    subjNum = s;
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end

    if contains(superorbasic,'Super')
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechan',superorbasic,'_LDAba_',trialtype,'Trials/']);
        if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
            mvpaoutputdir = fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']); % this is the direct MVPA output  
            group_mvpaoutput{s} = load(mvpaoutputdir).sl_map;
        end
    elseif contains(superorbasic,'Basic')
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechan',superorbasic,'1_LDAba_',trialtype,'Trials/']);
        save_dir2=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechan',superorbasic,'2_LDAba_',trialtype,'Trials/']);

         if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
            mvpaoutputdir = fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']); % this is the direct MVPA output  
            group_mvpaoutput{1,s} = load(mvpaoutputdir).sl_map;
         end
         if exist(fullfile([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
            mvpaoutputdir2 = fullfile([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']); % this is the direct MVPA output
            group_mvpaoutput{2,s} = load(mvpaoutputdir2).sl_map;
         end

    elseif contains(superorbasic,'Exemplar')
        time_radius = 0;        
        nchunks = 10;
        counter = 0;
            for p = 1:4
                pair = exemppairlist{p};
                for exempval1 = pair(1):pair(3) %one vs other schematic
                    for exempval2 = exempval1:pair(3) %one vs other schematic
                        if exempval1 ~= exempval2
                            counter = counter + 1;
                            save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/ExemplarDecoding1vs1_withactmap/',condition,'/noPCA_35Hz_100Hz_', num2str(nchunks), 'chunks_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_exemp',num2str(exempval1),num2str(exempval2),'_LDAba_',trialtype,'Trials/']);
                            if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
                            group_mvpaoutput{counter,s} = load([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat'],'sl_map').sl_map;
                            else
                                continue;
                            end
                        end
                    end
                end
            end
        end
    end

if contains(superorbasic,'Basic') || contains(superorbasic,'Exemplar')
    counter = 0;
    for s = 1:size(group_mvpaoutput,2)
        %remove empties
        subjgrp = group_mvpaoutput(:,s);
        subjgrp = subjgrp(~cellfun(@isempty, subjgrp));
        if ~isempty(subjgrp)
            counter = counter + 1;
            submatrix = cell2mat(cellfun(@(x) x.samples, subjgrp', 'UniformOutput', false)');
            avg{1,counter}.samples = mean(submatrix,1);
            avg{1,counter}.a = subjgrp{1}.a;
            avg{1,counter}.fa = subjgrp{1}.fa;
            avg{1,counter}.sa = subjgrp{1}.sa;
        end
    end
    n_subj = length(avg);
    group_mvpaoutput2 = avg;
else
    group_mvpaoutput2 = group_mvpaoutput(~cellfun(@isempty, group_mvpaoutput));
    n_subj = length(group_mvpaoutput2);
end

%% Group datasets together:
%some "non-unique" elements due to 0 being jittered, just use the
%parameters from subject 1 because they are mostly the same.
disp(n_subj);
ds_cell = cell(n_subj,1);
actual_dataset= cell(n_subj,1);
for p=1:n_subj
        actual_dataset{p,1}  = group_mvpaoutput2{p};
        actual_dataset{p,1}.sa.targets = 1;
        actual_dataset{p,1}.sa.chunks = p;
        %round values to the 3rd decimal because "0" is a very small number
        %that jitters around 0.
        actual_dataset{p,1}.a.fdim.values{1,1} = round(group_mvpaoutput2{p}.a.fdim.values{1,1},2);
end

[~,ds_group]=cosmo_mask_dim_intersect(actual_dataset);
ds_group = cosmo_stack(ds_group);

[~,maxind] = max(nanmean(ds_group.samples));
realpeaktime(c,aa) = timewindow(maxind);

%% Stats
allbalacc{c,aa} = ds_group.samples;

allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group, 'time',allow_clustering_over_time);

ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);
sig = find(abs(ds_z.samples)>1.96);
sigmask = abs(ds_z.samples)>1.96;

%% Plotting
plotcounter = plotcounter + 1;

if contains(analysistype,'Super')
    plotcolor = [50 130 246]./255;
elseif contains(analysistype,'Basic')
    plotcolor = [240 134 80]./255;
elseif contains(analysistype,'Exemplar')
    plotcolor = [58 224 213]./255;
end

clearvars temp
plottime = timewindow;
indstart = find(round(group_mvpaoutput2{1,1}.a.fdim.values{1,1},3) == plottime(1));
indend = find(round(group_mvpaoutput2{1,1}.a.fdim.values{1,1},3) == plottime(end));

for s = 1:n_subj
    temp(s,:) = group_mvpaoutput2{s}.samples;
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
alpha(0.25);

time_values = timewindow;

plotnum(aa) = plot(plottime,meansub,'Color',plotcolor,'LineWidth',1);

if aa == 1
    addval = 0.04;
elseif aa == 2
    addval = 0.02;
elseif aa == 3
    addval = 0.06;
end

if ~isempty(sig)
sigtime = time_values(sig);
onsettime(c,aa) = sigtime(1);
[~, peakind] = max(meansub);
peaktime_avg(c,aa) = time_values(peakind);

y_pos(1) = 0.49;
y_pos(2) = 0.48;%max(meansub)+addval+0.015;
y_pos(3) = 0.47;%max(meansub)+addval+0.015;

%make significant section bold:
if sum(sigmask) ~= 0
    boldvals = meansub; 
    boldvals(~sigmask) = NaN; 
    plot(plottime, boldvals, 'Color',plotcolor, 'LineWidth', 3.5); 
end
for i = 1:length(sig)
    line([sigtime(i)-0.005 sigtime(i)+0.005], [y_pos(aa) y_pos(aa)], ...
         'Color', plotcolor, 'LineWidth', 3);
end
end

xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
if contains(save_dir,'LDAauc')
    ylabel('AUC (chance=.5)');
elseif contains(save_dir,'LDAba')
    ylabel('Balanced Accuracy (chance=.5)');
else
    ylabel('Accuracy (chance=.5)');
end

xlabel('time (s)');
title(sprintf([condition, '-', trialtype, ' (N = ', num2str(n_subj),')']));
xlim([-0.4 0.8]);
ylim([0.46 0.545]);

%print first sig. time
if ~isempty(sig)
x_pos = sigtime(1);

firstsigtime_label(c,aa) = round(sigtime(1)*1000,0); 
text(x_pos, y_pos(aa)+0.006, [num2str(firstsigtime_label(c,aa)) 'ms'], 'FontSize', 10);
end


end

diff_acc = allbalacc{c,1} - allbalacc{c,2};
ds_group.samples = diff_acc;
allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group, 'time',allow_clustering_over_time);
ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0,'niter',1000);
sig_ind = find(abs(ds_z.samples)>1.96);
corrsigtimeall = plottime(sig_ind);
if ~isempty(corrsigtimeall)
for i = 1:length(corrsigtimeall)
line([corrsigtimeall(i)-0.005 corrsigtimeall(i)+0.005], [0.54 0.54], ...
         'Color', 'k', 'LineWidth', 3);
end

diffsig_label = round(corrsigtimeall(1)*1000,0); 
x_pos = plottime(sig_ind(1));
text(x_pos, 0.544, [num2str(diffsig_label) 'ms'], 'FontSize', 10);

end
end

%% Stats between categories
% column is N/M/L, row is super or basic
% allbalacc{1,1} vs allbalacc{1,2}
Supervsbasic{1} = allbalacc{1,1} - allbalacc{1,2};
% allbalacc{2,1} vs allbalacc{2,2}
Supervsbasic{2} = allbalacc{2,1} - allbalacc{2,2};
% allbalacc{3,1} vs allbalacc{3,2}
Supervsbasic{3} = allbalacc{3,1} - allbalacc{3,2};


%% Bootstrapped onsets:
bootstrapNUM = 20000;
additionalCOND = [];
condColor{1} = [50 130 246]./255;
condColor{2} = [240 134 80]./255;
condColor{3} = [58 224 213]./255;

clearvars bootstrapbycond_all
figure; set(gcf,'Color','w');
hold on
for c = 1:3
    condition = conditionList{c};
for aa = 1:length(superorbasicLIST)
    analysistype = superorbasicLIST{aa};
    temp = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/BootstrappingTwoTail',num2str(bootstrapNUM),'_LDA_',trialtype,'Trials_stats',additionalCOND,'/',condition,'/ds_z_',analysistype]),'ds_z');
    ds_z_all = temp.ds_z;
    for i = 1:length(ds_z_all)
        sigvals = find(abs(ds_z_all{i}.samples)>1.96);
        if ~isempty(sigvals)
            sig_onset(i) = timewindow(sigvals(1));
        else
            sig_onset(i) = nan;
        end
    end
    % sig_onset = sig_onset(~isnan(sig_onset));
    onsettime(c,aa) = nanmean(sig_onset);
    sig_onset_all{c,aa} = sig_onset;

    confidence_int(1,aa) = prctile(sig_onset, 2.5); %low
    confidence_int(2,aa) = prctile(sig_onset, 97.5); %high
  
end
CIbyCond_onset{c} = confidence_int; %CItdist; %confidence_int %yCI95

% yMeanbyCond_onset{c} = yMean;

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
title(conditionList{c});
% set(gca, 'XDir', 'reverse'); % Flips the x-axis
end

figure; set(gcf,'Color','w'); hold on
figcounter = 0;
for aa = 1:length(analysistypeLIST)
    disp(analysistypeLIST{aa});
for c = 1:3
    figcounter = figcounter + 1;
    subplot(3,3,figcounter)
    histogram(sig_onset_all{c,aa})
    xlim([-0.4 0.5]);
    fprintf('%.4f [%.4f %.4f]\n', onsettime(c,aa), CIbyCond_onset{c}(1,aa), CIbyCond_onset{c}(2,aa));
end
end

%% Peak latency test:
clearvars bootstrapbycond_all confidence_int
figure; set(gcf,'Color','w');
hold on
for c = 1:3
    condition = conditionList{c};
for aa = 1:length(analysistypeLIST)
    analysistype = analysistypeLIST{aa};
    clearvars bootstrap_data_all
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/BootstrappingTwoTail',num2str(bootstrapNUM),'_LDA_', ...
        trialtype,'Trials_stats',additionalCOND,'/',condition,'/bootstrap_balacc_data_',analysistype,'.mat']));
    for b = 1:bootstrapNUM
        [~,peakval(b)] = max(nanmean(bootstrap_data_all{b}.samples));
    end
    peaktimes = timewindow(peakval);
    peak_all{c,aa} = peaktimes;
    avgpeaktime(c,aa) = nanmean(peaktimes);
    avgpeaktime(c,aa) = peaktime_avg(c,aa);
    confidence_int(aa,:) = prctile(peaktimes, [2.5 97.5]); %low
end
CIbyCond_peak{c} = confidence_int; %CItdist; %confidence_int %yCI95

% yMeanbyCond_onset{c} = yMean;

subplot(1,3,c); hold on
bar(1,avgpeaktime(c,1),'FaceColor',condColor{1});
bar(2,avgpeaktime(c,2),'FaceColor',condColor{2});
bar(3,avgpeaktime(c,3),'FaceColor',condColor{3});
errorbar(1,avgpeaktime(c,1), avgpeaktime(c,1) - CIbyCond_peak{c}(1,1),CIbyCond_peak{c}(1,2) - avgpeaktime(c,1),'.','Color','k','Linewidth',1);
errorbar(2,avgpeaktime(c,2), avgpeaktime(c,2) - CIbyCond_peak{c}(2,1),CIbyCond_peak{c}(2,2) - avgpeaktime(c,2),'.','Color','k','Linewidth',1);
errorbar(3,avgpeaktime(c,3), avgpeaktime(c,3) - CIbyCond_peak{c}(3,1),CIbyCond_peak{c}(3,2) - avgpeaktime(c,3),'.','Color','k','Linewidth',1);
xlim([0 5]);
ylim([-0.4 0.8])
xticks(1:5);
xticklabels({'Super','Basic','Exemplar'});
xtickangle(90);
ylabel('Peak Latency (ms)');
title(conditionList{c});
% set(gca, 'XDir', 'reverse'); % Flips the x-axis
end

for aa = 1:length(analysistypeLIST)
    disp(analysistypeLIST{aa});
for c = 1:3
    fprintf('%.4f [%.4f %.4f]\n', avgpeaktime(c,aa), CIbyCond_peak{c}(aa,1), CIbyCond_peak{c}(aa,2));
end
end

%% Load permuted bootstrap data for onset/peak differences:
% conditions are scrambled.
% clc
% close all
conditionList = {'Normal','Masked','LSF'};
superorbasicLIST1 = {'Super','Super','Basic'};
superorbasicLIST2 = {'Basic','Exemplar','Exemplar'};

clearvars confidence_int superorbasic
figure; set(gcf,'Color','w');
figcounter = 0;
for c = 1:3
    for aa = 1:length(superorbasicLIST1)
    % Super - Basic
    % Super - Exemplar
    % Basic - Exemplar

    superorbasic{1} = superorbasicLIST1{aa};
    superorbasic{2} = superorbasicLIST2{aa};
    if contains(superorbasic{1},'Super') && contains(superorbasic{2},'Basic')
        aa1 = 1; aa2 = 2;
    elseif contains(superorbasic{1},'Super') && contains(superorbasic{2},'Exemplar')
        aa1 = 1; aa2 = 3;
    elseif contains(superorbasic{1},'Basic') && contains(superorbasic{2},'Exemplar')
        aa1 = 2; aa2 = 3;
    end

        
    condition = conditionList{c};
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/' ...
        'Null_BetweenCategory/BootstrappingTwoTail20000_LDA_',trialtype,['Trials_stats' ...
        '/'],condition,'/',superorbasic{1},'_',superorbasic{2},'/bootstrap_diff_onset',trialtype]),'bootstrap_diff_onset');
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/' ...
        'Null_BetweenCategory/BootstrappingTwoTail20000_LDA_',trialtype,['Trials_stats' ...
        '/'],condition,'/',superorbasic{1},'_',superorbasic{2},'/bootstrap_diff_peak',trialtype]),'bootstrap_diff_peak');

    bootstrap_diff_onset(isnan(bootstrap_diff_onset)) = 0; %set non-sig differences to 0
    bootstrap_diff_peak(isnan(bootstrap_diff_peak)) = 0; %set non-sig differences to 0

    confidence_int = prctile(bootstrap_diff_onset, [2.5 97.5]); %low
    CIbyCond_onset{c} = confidence_int;
    figcounter = figcounter + 1;
    realonsetdiff = firstsigtime_label(c,aa1)/1000 - firstsigtime_label(c,aa2)/1000;

    confidence_int = prctile(bootstrap_diff_peak, [2.5 97.5]); %low
    CIbyCond_peak{c} = confidence_int;

    subplot(3,3,figcounter)
    hold on
    histogram(bootstrap_diff_peak);
    realpeakdiff = realpeaktime(c,aa1) - realpeaktime(c,aa2);
    xline(realpeakdiff,'r');
    xline(confidence_int(1));
    xline(confidence_int(2));

    pval_onset(c,aa) = (sum(abs(bootstrap_diff_onset) >= abs(realonsetdiff)) + 1) / (length(bootstrap_diff_onset) + 1);
    pval_onset_onetail(c,aa) = (sum(bootstrap_diff_onset >= realonsetdiff) + 1) / (length(bootstrap_diff_onset) + 1);
    pval_peak(c,aa) = (sum(abs(bootstrap_diff_peak) >= abs(realpeakdiff)) + 1) / (length(bootstrap_diff_peak) + 1);
    end
end

for aa = 1:3
fprintf('Normal: Onset p-val = %.4f, Peak p-val = %.4f\n',  pval_onset(1,aa), pval_peak(1,aa));
fprintf('Masked: Onset p-val = %.4f, Peak p-val = %.4f\n',  pval_onset(2,aa), pval_peak(2,aa));
fprintf('LSF: Onset p-val = %.4f, Peak p-val = %.4f\n',  pval_onset(3,aa), pval_peak(3,aa));
end

%% Comparing confidence intervals to 0

% conditions are scrambled.
% clc
% close all
conditionList = {'Normal','Masked','LSF'};
superorbasicLIST1 = {'Super','Super','Basic'};
superorbasicLIST2 = {'Basic','Exemplar','Exemplar'};

clearvars confidence_int superorbasic
figure; set(gcf,'Color','w');
figcounter = 0;
for c = 1:3
    for aa = 1:length(superorbasicLIST1)
        
    % Super - Basic
    % Super - Exemplar
    % Basic - Exemplar

    superorbasic{1} = superorbasicLIST1{aa};
    superorbasic{2} = superorbasicLIST2{aa};
    if contains(superorbasic{1},'Super') && contains(superorbasic{2},'Basic')
        aa1 = 1; aa2 = 2;
    elseif contains(superorbasic{1},'Super') && contains(superorbasic{2},'Exemplar')
        aa1 = 1; aa2 = 3;
    elseif contains(superorbasic{1},'Basic') && contains(superorbasic{2},'Exemplar')
        aa1 = 2; aa2 = 3;
    end

        
    condition = conditionList{c};
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/' ...
        'DiffConfInt_BetweenCategory/BootstrappingTwoTail20000_LDA_',trialtype,['Trials_stats' ...
        '/'],condition,'/',superorbasic{1},'_',superorbasic{2},'/bootstrap_diff_onset',trialtype]),'bootstrap_diff_onset');
    load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/' ...
        'DiffConfInt_BetweenCategory/BootstrappingTwoTail20000_LDA_',trialtype,['Trials_stats' ...
        '/'],condition,'/',superorbasic{1},'_',superorbasic{2},'/bootstrap_diff_peak',trialtype]),'bootstrap_diff_peak');

    bootstrap_diff_onset(isnan(bootstrap_diff_onset)) = 0; %set non-sig differences to 0
    bootstrap_diff_peak(isnan(bootstrap_diff_peak)) = 0; %set non-sig differences to 0

    confidence_int = prctile(bootstrap_diff_onset, [2.5 97.5]); %low
    CIbyCond_onset{c} = confidence_int;

    confidence_int = prctile(bootstrap_diff_peak, [2.5 97.5]); %low
    figcounter = figcounter + 1;
    subplot(3,3,figcounter)
    hold on
    histogram(bootstrap_diff_peak);

    xline(confidence_int(1));
    xline(confidence_int(2));
    end
end

