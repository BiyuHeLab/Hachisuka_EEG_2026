%% For plotting CLMM results (outputs from R)
% This script will plot beta values as a function of time, and do
% cluster-size permuetation tests on the sig. timepoints.
% Note that the CLMM (cumulative-linked mix model) was fit for every 3 time bins (0.01 *3, or 30ms).
% Outcome was the PAS rating (ordinal variable), fixed effects were representational distances
% and conditions. Random effects were subjects.
% Last updated Feb 6 2026, Ayaka Hachisuka (ahachisu@gmail.com)

clear;
close all;
permanalysis = 'y'; % 'y' or 'n'
trialtype = 'ALL';
data = readtable(['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/Super_',trialtype,'_3timebin_signeddist.csv']);
dist_all{1} = data.beta1;
pval_real{1} = data.pval1;
t_real{1} = data.tval_dist;
data = readtable(['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/Basic_',trialtype,'_3timebin_signeddist.csv']);
dist_all{2} = data.beta1;
pval_real{2} = data.pval1;
t_real{2} = data.tval_dist;

categorylevel_list = {'Super','Basic'};

perm_file{1} = ['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/perms/Super/'];
perm_file{2} = ['/isilon/LFMI/VMdrive/Ayaka/EEG/Revision_round1/Fig5/AllCond/PASrating1/perms/Basic/'];

timeval = linspace(-0.4,0.8,240);
% Settings for permutatio nanalysis:
CDT = 0.05;
ClusterStats = 'ClusterSize'; % or 'SumPos'
n_perms = 10000;

if permanalysis == 'y'
    for s = 1:2
       
        % Step 1: Compute clusters in real data
        real_tvals = t_real{s};
        real_pvals = pval_real{s};
        clusters_orig = find_temporal_clusters(real_tvals, real_pvals, CDT);
        
        % Step 2: Compute cluster statistics for permutations
        maxClusterSizes = zeros(1, n_perms);
        maxStatSumPos   = zeros(1, n_perms);
        
        for i = 1:n_perms
            perm_data = readtable([perm_file{s},categorylevel_list{s},'_',trialtype,'_perm_',num2str(i),'.csv']);
            perm_tval = perm_data.tval_dist;
            perm_pval = perm_data.pval1;
            cluster_shuffled = find_temporal_clusters(perm_tval, perm_pval, CDT);
            maxClusterSizes(i) = cluster_shuffled.maxSize;
            maxStatSumPos(i)   = cluster_shuffled.maxStatSumPos;
        end
        
        % Step 3: Determine significance threshold
        if strcmp(ClusterStats, 'SumPos')
            maxStats = sort(maxStatSumPos, 'descend');
            CritVal = maxStats(round(0.05 * n_perms));  % Top 5% threshold
            SigClusters = find(clusters_orig.cluster_statSum > CritVal);
        elseif strcmp(ClusterStats, 'ClusterSize')
            maxStats = sort(maxClusterSizes, 'descend');
            % CritVal = maxStats(round(0.05 * n_perms));  % Top 5% threshold
            CritVal = prctile(maxStats, 95);
            SigClusters = find(clusters_orig.cluster_size > CritVal);
        end
        
        % Step 4: Get significant timepoints
        SigTimePoint = ismember(clusters_orig.cluster_timecourse, SigClusters);
        corrsigtimeall{s} = timeval(SigTimePoint);
    end
else
    ntimepoints = 240;
    alpha_corrected = CDT / ntimepoints;
    alpha = 0.05;
    for s = 1:4
        orig_pval = pval_real{s};
        % sig_ind = orig_pval < alpha; % < alpha_corrected;
        pvals_fdr = mafdr(orig_pval, 'BHFDR', true);
        sig_ind = pvals_fdr < alpha;
        corrsigtimeall{s} = timeval(sig_ind);
    end
end

%% Plotting
plotcolor{1} = [50 130 246]./255;
plotcolor{2} = [240 134 80]./255;
ypos(1) = 0.04; ypos(2) = 0.045;
figure; set(gcf,'Color','w'); hold on
for i = 1:2
corrsigtime = corrsigtimeall{i};
sigmask = ismember(round(timeval, 4), round(corrsigtime, 4));
plot(timeval,dist_all{i},'Color',plotcolor{i},'LineWidth',1);

%make significant section bold:
if sum(sigmask) ~= 0
    boldvals = dist_all{i}; 
    boldvals(~sigmask) = NaN; 
    plot(timeval, boldvals, 'Color',plotcolor{i}, 'LineWidth', 3.5); 
end
grid on
ylabel('Beta (distance)');
xlabel('Time (s)');
xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0, '--', 'LineWidth', 1, 'Color', 'k');
% for t = 1:length(corrsigtime)
%     line([corrsigtime(t)-0.005 corrsigtime(t)+0.005], [max(dist_all{i})+0.005 max(dist_all{i})+0.005], ...
%          'Color', plotcolor{i}, 'LineWidth', 3);
% end
sig_idx = find(sigmask);% Find clusters of consecutive significant points
d_sig = diff([0 sigmask 0]);  
cluster_starts = find(d_sig == 1);
cluster_ends = find(d_sig == -1) - 1;
y_pos1 = [0.07 0.05 0.05 0.05];
y_pos2 = [-0.04 -0.05 -0.05 -0.05];
hold on;
data = dist_all{i};
% for c = 1:length(cluster_starts)
%     idx_range = cluster_starts(c):cluster_ends(c);
%     cluster_data = data(idx_range);
%     cluster_time = timeval(idx_range);
% 
%     % Decide cluster polarity by mean data sign
%     if mean(cluster_data) >= 0 
%         % Positive cluster: bar above max data value
%         y_bar = y_pos1(i); %max(cluster_data) + 0.03*range(data);
%     else
%         % Negative cluster: bar below min data value
%         y_bar = y_pos2(i); %min(cluster_data) - 0.03*range(data);
%     end
% 
%     % Plot horizontal line over entire cluster time range
%     plot([cluster_time(1), cluster_time(end)], [y_bar y_bar], 'Color', plotcolor{i}, 'LineWidth', 3);
% end

end
% ylim([-0.4 0.5]);
xlim([-0.4 0.8]);
title('Predicting categorization accuracy (PAS = 1 trials only)');
