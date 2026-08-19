%% Bootstrapping 1000 times.

clear;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa');
%% Load each subject:
numSubjects = 33;
subjectslist = setdiff(1:numSubjects,[22 32 33]);
bootstrapNUM = 20000;
conditionLIST = {'Normal','Masked','LSF'};
exemppairlist{1} = 1:3;
exemppairlist{2} = 4:6;
exemppairlist{3} = 7:9;
exemppairlist{4} = 10:12;

numCores = 18; %feature('numcores')-2;

for cc = 1:length(conditionLIST)
condition = conditionLIST{cc};
superorbasicLIST = {'Super','Basic','Exemplar'};
analysistypeLIST = {'Super','Basic','Exemplar'};
% superorbasicLIST = {'Exemplar'};
% analysistypeLIST = {'Exemplar'};
trialtype = 'ALL';
time_radius = 0;        
trial_bin_num = 4;
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
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/',condition,'/LDA/noPCA_35Hz_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_',num2str(numrepeats),'repeats_morechan',superorbasic,'_LDAba_',trialtype,'Trials/']);
        if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
            mvpaoutputdir = fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']); % this is the direct MVPA output  
            group_mvpaoutput{sub_counter} = load(mvpaoutputdir).sl_map;
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
        trial_bin_num = 4;
        nchunks = 10;
        counter = 0;
         for p = 1:4
            pair = exemppairlist{p};
            for exempval1 = pair(1):pair(3) %one vs other schematic
                    for exempval2 = exempval1:pair(3) %one vs other schematic
                        if exempval1 ~= exempval2                    counter = counter + 1;
                        subjNum = s;
                        if length(num2str(subjNum)) == 1 || ...
                                (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
                            subjNumStr = ['0' num2str(subjNum)];
                        else
                            subjNumStr = num2str(subjNum);
                        end
                        
                        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/ExemplarDecoding1vs1_withactmap/',condition,'/noPCA_35Hz_100Hz_', num2str(nchunks), 'chunks_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_exemp',num2str(exempval1),num2str(exempval2),'_LDAba_',trialtype,'Trials/']);
                        if exist(fullfile([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat']),'file')
                        group_mvpaoutput{counter,s} = load([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput.mat'],'sl_map').sl_map;
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
            avg{1,counter}.samples = nanmean(submatrix,1);
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
clearvars ds_group actual_dataset
ds_cell = cell(n_subj,1);
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
% n_subj = size(ds_group.samples,1);
%% Bootstrapping

% timewindow = linspace(-0.4,0.8,240);
timewindow = -0.4:0.005:0.8;
clearvars temp
plottime = timewindow;
indstart = find(round(group_mvpaoutput2{1,1}.a.fdim.values{1,1},2) == plottime(1));
indend = find(round(group_mvpaoutput2{1,1}.a.fdim.values{1,1},2) == plottime(end));
indstart = indstart(1); %to make sure it's 1 value
indend = indend(end); %to make sure it's 1 value

for s = 1:n_subj
    decodingacc_bysub(s,:) = group_mvpaoutput2{s}.samples(indstart:indend);
end

%bootstrapping for statistics:
allow_clustering_over_time = true;

bootstrap_sig_onset = zeros(1,bootstrapNUM);
clearvars bootstrapgroup
parfor (b = 1:bootstrapNUM, numCores)
% for b = 1:bootstrapNUM
    bootstrap_indices = randi(n_subj, [1, n_subj]); %sample from 29 subjects, with replacement
    bootstrapacc = decodingacc_bysub(bootstrap_indices,:);

    bootstrapgroup = [];
    bootstrapgroup.samples = bootstrapacc;
    bootstrapgroup.sa = ds_group.sa;
    bootstrapgroup.fa = ds_group.fa;
    bootstrapgroup.a = ds_group.a;
    bootstrapgroup.a.fdim.values{1,1} = ds_group.a.fdim.values{1,1}(indstart:indend);
    bootstrapgroup.fa.center_ids = ds_group.fa.center_ids(indstart:indend);
    bootstrapgroup.fa.time = ds_group.fa.time(indstart:indend)
    
    nbrhood = cosmo_cluster_neighborhood(bootstrapgroup, 'time',allow_clustering_over_time);
    % ds_z{b} = cosmo_montecarlo_cluster_stat(bootstrapgroup,nbrhood,'h0_mean',0.5,'niter',1000,'progress',[]);
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

savedir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/BootstrappingTwoTail',num2str(bootstrapNUM),'_LDA_',trialtype,'Trials_stats/',condition,'/'];
if ~exist(savedir,'file'); mkdir(savedir); end
save(fullfile([savedir,'bootstrap_sigval_',analysistype]),'bootstrapbycond');
save(fullfile([savedir,'bootstrap_balacc_data_',analysistype]),'bootstrap_data_all');
save(fullfile([savedir,'ds_z_',analysistype]),'ds_z');

end
end

%% Plot onsets
% clearvars plotcolor
% 
% superorbasicLIST = {'Exemplar','Basic','Super'};
% analysistypeLIST = {'Exemplar','Basic','Super'};
% superorbasicLIST = {'Super','Basic'};
% analysistypeLIST = {'Super','Basic'};
% conditionList = {'Normal','Masked','LSF'};
% 
% plotcolor(1,:) = [255 73 166]./255;
% plotcolor(2,:) = [135 97 255]./255;
% plotcolor(3,:) = [58 224 213]./255;
% 
% 
% figure; set(gcf,'Color','w');
% for c = 1:length(conditionList)
%     condition = conditionList{c};
% for aa = 1:length(superorbasicLIST)
%     analysistype = analysistypeLIST{aa};
%     temp = load(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/Bootstrapping',num2str(bootstrapNUM),'_LDA_ALLTrials_stats/',condition,'/bootstrap_sigval_',analysistype]),'bootstrapbycond');
%     bootstrapbycond_all(aa,:) = temp.bootstrapbycond;
% end
%     yMean(1) = nanmean(bootstrapbycond_all(1,:));
%     yMean(2) = nanmean(bootstrapbycond_all(2,:));
%     yMean(3) = nanmean(bootstrapbycond_all(3,:));
%     yMeanbyCond{c} = yMean;
%     N = size(bootstrapbycond_all,2);
%     ySEM(1) = nanstd(bootstrapbycond_all(1,:))/sqrt(N);    
%     ySEM(2) = nanstd(bootstrapbycond_all(2,:))/sqrt(N);    
%     ySEM(3) = nanstd(bootstrapbycond_all(3,:))/sqrt(N);    
% 
%     CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
%     yCI95(:,1) = bsxfun(@times, ySEM(1), CI95(:));
%     yCI95(:,2) = bsxfun(@times, ySEM(2), CI95(:));
%     yCI95(:,3) = bsxfun(@times, ySEM(3), CI95(:));
%     CIbyCond{c} = round(yCI95+yMean,2);
% end
% 
% 
% 
% subplot(1,3,c)
% hold on
% for cat = 1:3
% h(cat) = boxchart(repmat(cat,1,size(bootstrapbycond_all,2)),bootstrapbycond_all(cat,:),'MarkerStyle','none');
% h = plot(yMean+yCI95);
% 
% if cat == 3
%     plotcolor = [0 18 154]./255;
% elseif cat == 2
%     plotcolor = [115 43 245]./255;
% elseif cat == 1
%     plotcolor = [58 224 213]./255;
% end
% 
% h(cat).BoxFaceColor = plotcolor;
% ylabel('Decoding onset (ms)');
% title(condition);
% ylim([-100 600]);
% end
% xticks = 1:3;
% xticklabels(superorbasicLIST);
% 
% 
% mean1 = nanmean(bootstrapbycond(1,:));
% mean2 = nanmean(bootstrapbycond(2,:));
% [H(1), pValue(1), CI(1,:), stats{1}] = ttest2(bootstrapbycond_all(1,:), bootstrapbycond_all(2,:), 'Vartype','unequal');
% [H(2), pValue(2), CI(2,:), stats{2}] = ttest2(bootstrapbycond_all(1,:), bootstrapbycond_all(3,:), 'Vartype','unequal');
% [H(3), pValue(3), CI(3,:), stats{3}] = ttest2(bootstrapbycond_all(2,:), bootstrapbycond_all(3,:), 'Vartype','unequal');
% 
% [~, p_value_ttest(c), ~, stats] = ttest2(bootstrapbycond_all(1,:), bootstrapbycond_all(2,:));
% p_value_wilcoxon(c) = ranksum(bootstrapbycond_all(1,:), bootstrapbycond_all(2,:));
% 
% [~, p_value_ttest2, ~, stats] = ttest2(bootstrapbycond_all(1,:), bootstrapbycond_all(3,:));
% p_value_wilcoxon2 = ranksum(bootstrapbycond_all(1,:), bootstrapbycond_all(3,:));
% 
% [~, p_value_ttest3, ~, stats] = ttest2(bootstrapbycond_all(2,:), bootstrapbycond_all(3,:));
% p_value_wilcoxon3 = ranksum(bootstrapbycond_all(2,:), bootstrapbycond_all(3,:));
% 
% 
% function consec_onset = find_onset(zvalues)
%     consec_count = 0;
%     consec_onset = NaN;
%     for t = 1:length(ds_z.samples)
%         if ds_z.samples > 1.64
%             consec_count = consec_count + 1;
%         else
%             consec_count = 0;
%         end
% 
%         if consec_count == 10 %once minimum is reached
%             consec_onset = t - 9;
%             break;
%         end
%     end
% end
