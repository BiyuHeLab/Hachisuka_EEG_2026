%% Calculate & Plot Parallelism Score from CCGP

% Last updated June 20 2025, Ayaka Hachisuka

clear;
close all
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');

%% Choose analysis type:
analysistypeLIST = {'supervsbasicall'};
yvals = [0.12 0.115 0.1];
%% Initialize parameters:
numSubjects = 33;
numrepeats = 100;
subjectslist = setdiff(1:numSubjects,[22 32 33]);
cond_all = {'Normal','Masked','LSF'};
pairlist1 = {'1A','2A'};
pairlist2 = {'1B','2B'};
timevals = -0.4:0.005:0.8;

for a = 1:size(analysistypeLIST,2)
analysistype = analysistypeLIST{a};

basiccat_list = [];
totalpair = 2;
if strcmp(analysistype,'supervsbasicall')
    analysistype2 = '_BasicABCDALL';
    totalpair = 2;
elseif strcmp(analysistype,'supervsbasicdetect')
    analysistype2 = '_BasicABCDDetect';
    totalpair = 2;
elseif strcmp(analysistype,'basiccontrol')
    analysistype2 = '_BasicControlALL';
    totalpair = 1;
end

startIND = find(round(timevals,3) == 0);
endIND = find(round(timevals,3) ==0.3);

pairlist1 = [1 1 2 2];
pairlist2 = [1 1 2 2];
comblist1a = [1 2 1 2];
comblist1b = [2 1 2 1];
comblist2a = [2 1 2 1];
comblist2b = [1 2 1 2];
%% Calculate PS
for c = 1:3

    condition = cond_all{c};
    savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/resultswith200Hzdata/',condition,analysistype2]);
    %savedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CCGP_inputs/resultswith200HzdataV2/',condition, analysistype2, '_balance_train_only']);

    subcounter = 0;
    for s = subjectslist
        subcounter = subcounter + 1;
        paircounter = 0;
        % temp = nan(4,241);
        for ii = 1:length(pairlist1)
            pairnum1 = pairlist1(ii);
            pairnum2 = pairlist2(ii);
            comb1a = comblist1a(ii);
            comb1b = comblist1b(ii);
            comb2a = comblist2a(ii);
            comb2b = comblist2b(ii);
            if exist(fullfile([savedir,'/sub', num2str(s), '_pair',num2str(pairnum1),'_combo_',num2str(comb1a),num2str(comb1b),'_',num2str(numrepeats),'repeats_ds_searchlight_result.mat']))...
                && exist(fullfile([savedir,'/sub', num2str(s), '_pair',num2str(pairnum2),'_combo_',num2str(comb2a),num2str(comb2b),'_',num2str(numrepeats),'repeats_ds_searchlight_result.mat']))
                    paircounter = paircounter + 1;
                    load(fullfile([savedir,'/sub', num2str(s), '_pair',num2str(pairnum1),'_combo_',num2str(comb1a),num2str(comb1b),'_',num2str(numrepeats),'repeats_ds_searchlight_result']),'ds_searchlight_result');
                    hyperplane1 = ds_searchlight_result.hyperplane;
                    load(fullfile([savedir,'/sub', num2str(s), '_pair',num2str(pairnum2),'_combo_',num2str(comb2a),num2str(comb2b),'_',num2str(numrepeats),'repeats_ds_searchlight_result']),'ds_searchlight_result');
                    hyperplane2 = ds_searchlight_result.hyperplane;
                    for t = 1:length(timevals)
                        a = hyperplane1(:,t);
                        b = hyperplane2(:,t);
                        hyper_temp(paircounter,t) = (dot(a,b)/(norm(a)*norm(b))); % how you calculate cosine similarlity. the lower the number, the more parallel
                        % alternative (same result): v = [a';b']; hyper_temp(paircounter,t) = 1-pdist(v,'cosine');
                    end
            end
        end
        if totalpair > 1
            ps_temp1 = mean([hyper_temp(1,:);hyper_temp(2,:)]);
            ps_temp2 = mean([hyper_temp(3,:);hyper_temp(4,:)]);
            ps_temp = mean([ps_temp1;ps_temp2]);
            baseline = ps_temp(1:81); %ind = 81 is time = 0;
            mean_bas = mean(baseline);
            average_ps_baselinecorr = ps_temp - mean_bas;

            ps{c}(subcounter,:) = average_ps_baselinecorr; %smoothdata(average_ps_baselinecorr,'sgolay',20);
            ps_sem{c}(subcounter,:) = nanstd(average_ps_baselinecorr)/sqrt(subcounter);
        else
            ps_temp = mean(hyper_temp);
            baseline = ps_temp(1:81); %ind = 81 is time = 0;
            mean_bas = mean(baseline);
            average_ps_baselinecorr = ps_temp - mean_bas;

            ps{c}(subcounter,:) = average_ps_baselinecorr; %smoothdata(average_ps_baselinecorr,'sgolay',20);
            ps_sem{c}(subcounter,:) = nanstd(average_ps_baselinecorr)/sqrt(subcounter);
        end
    end
    average_ps(c,:) = nanmean(ps{c});
    % baseline_stats(c) = nanmean(nanmean(ps{c}(:,1:81)));
    sem_ps(c,:) = nanstd(ps_sem{c})/sqrt(subcounter); %nanstd(ps_sem{c})/sqrt(subcounter);
    end
end

%% save for BF stats:
% save ps{c}

%% Statistics
for c = 1:3
    for p = 1:subcounter
        ps_all{p,1}.samples = ps{c}(p,:);
        ps_all{p,1}.sa.targets = 1;
        ps_all{p,1}.sa.chunks = p;
        ps_all{p,1}.sa.labels = {'parallelism_score'};
        ps_all{p,1}.a.fdim.values{1,1} = linspace(-0.4,0.8,241);
        ps_all{p,1}.a.fdim.labels = {'time'};
        ps_all{p,1}.a.meeg.samples_field = 'trial';
        ps_all{p,1}.fa.center_ids = 1:241;
        ps_all{p,1}.fa.time = 1:241;
    end
    [~,ds_group]=cosmo_mask_dim_intersect(ps_all);
    ds_group = cosmo_stack(ds_group);
    allow_clustering_over_time = true;
    nbrhood = cosmo_cluster_neighborhood(ds_group, 'time',allow_clustering_over_time);
    % ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0,'niter',10000);
    ds_z = cosmo_montecarlo_cluster_stat(ds_group,nbrhood,'h0_mean',0,'niter',1000,'cluster_stat','maxsize','p_uncorrected',0.05);
    sig = find(abs(ds_z.samples)>1.64);
    if ~isempty(sig)
        corrsigtimeall{c} = timevals(sig);
    else
        corrsigtimeall{c} = [];
    end
end

% Rearrange dat afor visualization.
for c = 1:3
    temp_ps = average_ps(c,:);
    average_ps_baselinecorr(c,:) = smoothdata(temp_ps,'sgolay',10);
end

%% Plotting
condColors = {[0.9290, 0.6940, 0.1250],[0.6350, 0.0780, 0.1840],[0.4660 0.6740 0.1880]};
timewindow = linspace(-0.4,0.8,241);

figure; set(gcf,'Color','w');
hold on
for p = 1:3
    % subplot(3,1,p); hold on
    plot(timevals,average_ps_baselinecorr(p,:),'Color',condColors{p},'LineWidth',1.5);
    x = timewindow;
    curve1 = average_ps_baselinecorr(p,:) + sem_ps(p,:);
    curve2 = average_ps_baselinecorr(p,:) - sem_ps(p,:);
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    h = fill(x2, inBetween,condColors{p});
    set(h, 'EdgeColor', 'none'); % Remove black border
    alpha(0.25);
    
for c = 1:3
    sigtime = corrsigtimeall{c};
    if ~isempty(sigtime)
        for i = 1:length(sigtime)
            % Draw a horizontal line across the full x-axis at each significant value
            line([sigtime(i)-0.005 sigtime(i)+0.005], [yvals(c) yvals(c)], ...
                 'Color', condColors{c}, 'LineWidth', 3);
        end
    end
end
grid 'on'
xline(0, '--', 'LineWidth', 1, 'Color', 'k');
yline(0, '--', 'LineWidth', 1, 'Color', 'k');
xlabel('Time(s)');
ylabel('Parallelism score');
xlim([-0.4 0.9]);
ylim([-0.065 0.15]);
set(gca, 'FontName', 'Arial', 'FontSize', 25)
set(gca, 'YGrid', 'on');
end

