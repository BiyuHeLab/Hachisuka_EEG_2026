%% Group analysis for undetected, all conditions:

% Last updated Feb 6 2026, Ayaka Hachisuka

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
%% Load each subject:
numSubjects = 33;
% conditionList = {'Normal','Masked','LSF'};
plotcounter = 0;

superorbasicLIST = {'Super','Basic'};
analysistypeLIST = {'Super','Basic'};

trialtype = 'Undetect';
time_radius = 0;        
trial_bin_num = 1;
numrepeats = 300;
OCCIP = 0;
run_pca = 0;
subjectslist = setdiff(1:numSubjects,[22 32 33]);

exemppairlist{1} = 1:3;
exemppairlist{2} = 4:6;
exemppairlist{3} = 7:9;
exemppairlist{4} = 10:12;
condition = 'AllCond';
% subplot(3,1,c);
% hold on
figure; set(gcf,'Color','w');
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
       save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/CombinedConditions/AllCond/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechanSuper_LDAba_',trialtype,'Trials/']);
        if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
            mvpaoutputdir = fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']); % this is the direct MVPA output  
            group_mvpaoutput{s} = load(mvpaoutputdir).sl_map;
        else
            continue;
        end
    elseif contains(superorbasic,'Basic')
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/CombinedConditions/AllCond/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechanBasic1_LDAba_',trialtype,'Trials/']);
        save_dir2=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/CombinedConditions/AllCond/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechanBasic2_LDAba_',trialtype,'Trials/']);

         if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
            mvpaoutputdir = fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']); % this is the direct MVPA output  
            group_mvpaoutput{1,s} = load(mvpaoutputdir).sl_map;
         end
         if exist(fullfile([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
            mvpaoutputdir2 = fullfile([save_dir2,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']); % this is the direct MVPA output
            group_mvpaoutput{2,s} = load(mvpaoutputdir2).sl_map;
         else
            continue;
         end

    elseif contains(superorbasic,'Exemplar')
        time_radius = 0;        
        counter = 0;
            for p = 1:4
                pair = exemppairlist{p};
                for exempval1 = pair(1):pair(3) %one vs other schematic
                    for exempval2 = exempval1:pair(3) %one vs other schematic
                        if exempval1 ~= exempval2
                            counter = counter + 1;
                            save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/ExemplarDecoding1vs1/CombinedCond/noPCA_35Hz_100Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_exemp',num2str(exempval1),num2str(exempval2),'_LDAba_',trialtype,'Trials/']);
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
            grp{1,counter}.samples = mean(submatrix,1);
            grp{1,counter}.a = subjgrp{1}.a;
            grp{1,counter}.fa = subjgrp{1}.fa;
            grp{1,counter}.sa = subjgrp{1}.sa;
        end
    end
    n_subj = length(grp);
    group_mvpaoutput2 = grp;
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

%% Stats

allow_clustering_over_time = true;
nbrhood = cosmo_cluster_neighborhood(ds_group, 'time',allow_clustering_over_time);

ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0.5,'niter',10000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);
sig = find(ds_z.samples>1.96);
sigmask = abs(ds_z.samples)>1.96;

%% Plotting
plotcounter = plotcounter + 1;

if contains(analysistype,'Super')
    plotcolor = [50 130 246]./255;
elseif contains(analysistype,'Basic')
    plotcolor = [240 134 80]./255;
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
title(sprintf(['AllCond-', trialtype, ' (N = ', num2str(n_subj),')']));
xlim([-0.4 0.8]);
ylim([0.46 0.56]);

%print first sig. time
if ~isempty(sig)
x_pos = sigtime(1);

firstsigtime_label(c,aa) = round(sigtime(1)*1000,0); 
text(x_pos, y_pos(aa)+0.006, [num2str(firstsigtime_label(c,aa)) 'ms'], 'FontSize', 10);
end

end
