%% Master Script for OSOM analysis (Experiments 1 & 2)
% Plots Super and Basic Category performance for each condition.
% Plots PAS ratings for each condition.
% Script adapted from Thomas Baumgarten
% Last updated June 20 2025, Ayaka Hachisuka (ahachisu@gmail.com)

clear;
load('/isilon/LFMI/VMdrive/Ayaka/EEG/OSOM_fromThomas/Group_n16_Behav_DescrStats.mat');

% Note, trials included 435 target image trials and 141 catch trials
% pas_target & pas_catch are formatted as: (trial x sub x pas x condition 1-3)

%% CHOOSE Exp1 or Exp2:
experiment = 'Exp1'; % 'Exp1' or 'Exp2'

if strcmp(experiment,'Exp1')
    subjectslist = 1:16;
    numSubjects = 16;
    cond1 = 'Normal';
    cond2 = 'Masked';
    cond3 = 'MaskedMoved';
    cond1name = 'Normal';
    cond2name = 'OSM';
    cond3name = 'OSOM';
    condColors = {[0.9290, 0.6940, 0.1250],[0.8500, 0.3250, 0.0980],[0.6350, 0.0780, 0.1840]};
    foldername = 'OSOM_Exp1';
elseif strcmp(experiment,'Exp2')
    subjectslist = 17:25;
    numSubjects = 9;
    cond1 = 'OSM';
    cond2 = 'OSOMstable';
    cond3 = 'OSOMmoved';
    cond1name = cond1;
    cond2name = cond2;
    cond3name = cond3;
    condColors = {[0.8500, 0.3250, 0.0980],[0.6350, 0.0780, 0.1840],[112/255, 90/255, 68/255]};
    foldername = 'OSOM_Exp2';
end

%% Initialize variables
cond_list = {cond1,cond2,cond3};
pas_target = cell(numSubjects,1);
pas_catch = cell(numSubjects,1);

%% MEAN PAS CALCULATION:
levelist = {'Super','Basic'};
pas_all = cell(3,1); %[];
pas_supercorr = cell(3,1); %[];
pas_basiccorr = cell(3,1); %[];
for cond = 1:3
    pas_all{cond} = cell(numSubjects,1);
    pas_supercorr{cond} = cell(numSubjects,1);
    pas_basiccorr{cond} = cell(numSubjects,1);
end

subcounter = 0;
for i_sub = subjectslist
    subcounter = subcounter + 1;
    pasresp = table2array(DescrStats{1,i_sub}.Trial(:,13));
    condition = table2array(DescrStats{1,i_sub}.Trial(:,5));
    targetvscatch = table2array(DescrStats{1,i_sub}.Trial(:,4));

    supercorrect = table2array(DescrStats{1,i_sub}.Trial(:,14));
    basiccorrect = table2array(DescrStats{1,i_sub}.Trial(:,15));
    basiccorrect(basiccorrect<0)=NaN;

    pas_all{1}{subcounter} = pasresp(contains(condition,cond1) & contains(targetvscatch,'Target'));
    pas_supercorr{1}{subcounter} = pasresp(supercorrect==1 & contains(condition,cond1) & contains(targetvscatch,'Target'));
    pas_basiccorr{1}{subcounter} = pasresp(basiccorrect==1 & contains(condition,cond1) & contains(targetvscatch,'Target'));
    pas_all{2}{subcounter} = pasresp(contains(condition,cond2) & contains(targetvscatch,'Target'));
    pas_supercorr{2}{subcounter} = pasresp(supercorrect==1 & contains(condition,cond2) & contains(targetvscatch,'Target'));
    pas_basiccorr{2}{subcounter} = pasresp(basiccorrect==1 & contains(condition,cond2) & contains(targetvscatch,'Target'));
    pas_all{3}{subcounter} = pasresp(contains(condition,cond3) & contains(targetvscatch,'Target'));
    pas_supercorr{3}{subcounter} = pasresp(supercorrect==1 & contains(condition,cond3) & contains(targetvscatch,'Target'));
    pas_basiccorr{3}{subcounter} = pasresp(basiccorrect==1 & contains(condition,cond3) & contains(targetvscatch,'Target'));
end

for cond = 1:3
    for i_sub = 1:numSubjects
        pas_allsub(i_sub,cond) = nanmean(pas_all{cond}{i_sub});
        pas_allsub_supercorr(i_sub,cond) = nanmean(pas_supercorr{cond}{i_sub});
        pas_allsub_basiccorr(i_sub,cond) = nanmean(pas_basiccorr{cond}{i_sub});
    end
end

%saving for stats:
clearvars pasbytrial2 condbytrial2 sub2 pasbytrial condbytrial sub stats_table
for i_sub = 1:numSubjects
    for cond = 1:3
    temptriallist = pas_all{cond}{i_sub};
    pasbytrial2{cond} = temptriallist;
    condbytrial2{cond} = ones(length(temptriallist),1)*cond;
    sub2{cond} = repmat(i_sub,length(temptriallist),1);
    end
    pasbytrial{i_sub} = cell2mat(pasbytrial2(:));
    condbytrial{i_sub} = cell2mat(condbytrial2(:));
    sub{i_sub} = cell2mat(sub2(:));
end

%% Plot Mean PAS
analysistypelist = {'PAS with all trials','PAS with correct super trials','PAS with correct basic trials'};
figure; set(gcf,'Color','w');
hold on
for aa = 1:length(analysistypelist)
    analysistype = analysistypelist{aa};
if contains(analysistype,'PAS with all trials')
    normal_pas = pas_allsub(:,1);
    masked_pas = pas_allsub(:,2);
    lsf_pas = pas_allsub(:,3);
    % save('/isilon/LFMI/VMdrive/Ayaka/EEG/Behavior_stats/pas_allsub_alltrials.mat','pas_allsub');
elseif contains(analysistype,'PAS with correct super trials')
    normal_pas = pas_allsub_supercorr(:,1);
    masked_pas = pas_allsub_supercorr(:,2);
    lsf_pas = pas_allsub_supercorr(:,3);
    % save('/isilon/LFMI/VMdrive/Ayaka/EEG/Behavior_stats/pas_allsub_supercorrtrials.mat','pas_allsub_supercorr');
elseif contains(analysistype,'PAS with correct basic trials')
    normal_pas = pas_allsub_basiccorr(:,1);
    masked_pas = pas_allsub_basiccorr(:,2);
    lsf_pas = pas_allsub_basiccorr(:,3);
    % save('/isilon/LFMI/VMdrive/Ayaka/EEG/Behavior_stats/pas_allsub_basiccorrtrials.mat','pas_allsub_basiccorr');
end

all_pas = [normal_pas masked_pas lsf_pas];

subplot(3,1,aa)
hold on
offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows                 
for r = 1:3              
    xPos = repmat(offsets(r),1,numSubjects);
    boxchart(xPos, all_pas(:,r),'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off');  
end

xticklabels([{''},{''},{''},cond1name,cond2name,cond3name,{''}]);
ylabel('mean PAS');
ylim([0 2.2]);
% xlabel('Condition');
set(gca, 'YGrid', 'on');

clearvars xval
               
for r = 1:3               
    xPos = repmat(offsets(r),1,numSubjects);
    individDataPts = all_pas(:,r);
    s = swarmchart(xPos,individDataPts,30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', condColors{r}); %,'*','LineWidth',2,'Color','k');
    individDataPts2(r,:) = individDataPts;
    s.XJitter = 'rand';
    s.XJitterWidth = 0.1;
    s.MarkerFaceAlpha = .7;
    xval(r,:) = s.XData;
end


for i_sub = 1:numSubjects
    plot([xval(1,i_sub),xval(2,i_sub),xval(3,i_sub)],[individDataPts2(1,i_sub), individDataPts2(2,i_sub), individDataPts2(3,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
end

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows                 
for r = 1:3              
    xPos = repmat(offsets(r),1,numSubjects);
    boxchart(xPos, all_pas(:,r),'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off');  
end

title(analysistype);
end

%% Superordinate & Basic Category performance (Target trials only)

clearvars pasresp condition targetvscatch supercorrect basiccorrect
scounter = 0;
for i_sub = subjectslist
    scounter = scounter + 1;
    pasresp = table2array(DescrStats{1,i_sub}.Trial(:,13));
    condition = table2array(DescrStats{1,i_sub}.Trial(:,5));
    targetvscatch = table2array(DescrStats{1,i_sub}.Trial(:,4));

    supercorrect = table2array(DescrStats{1,i_sub}.Trial(:,14));
    basiccorrect = table2array(DescrStats{1,i_sub}.Trial(:,15));
    basiccorrect(basiccorrect<0)=NaN; % if super is incorrect, basiccorrect is recorded as -1. Make this 0 if you want to count them as incorrect. To ignore, keep it NaN.
    for c = 1:3
    condtype = cond_list{c};
    supercorrect_all{c}{scounter} = supercorrect(contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(supercorrect));
    supercorrect_undetect{c}{scounter} = supercorrect(pasresp == 0 & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(supercorrect));
    supercorrect_detect{c}{scounter} = supercorrect((pasresp == 1 | pasresp == 2 | pasresp == 3) & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(supercorrect));
    supercorrect_unrec{c}{scounter} = supercorrect((pasresp == 0 | pasresp == 1) & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(supercorrect));
    supercorrect_rec{c}{scounter} = supercorrect((pasresp == 2 | pasresp == 3) & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(supercorrect));

    basiccorrect_all{c}{scounter} = basiccorrect(contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(basiccorrect));
    basiccorrect_supercorr{c}{scounter} = basiccorrect(supercorrect == 1 & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(basiccorrect));
    basiccorrect_undetect{c}{scounter} = basiccorrect(supercorrect == 1 & pasresp == 0 & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(basiccorrect));
    basiccorrect_detect{c}{scounter} = basiccorrect(supercorrect == 1 & (pasresp == 1 | pasresp == 2 | pasresp == 3) & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(basiccorrect));
    basiccorrect_unrec{c}{scounter} = basiccorrect(supercorrect == 1 & (pasresp == 0 | pasresp == 1) & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(basiccorrect));
    basiccorrect_rec{c}{scounter} = basiccorrect(supercorrect == 1 & (pasresp == 2 | pasresp == 3) & contains(condition,condtype) & contains(targetvscatch,'Target'));% & ~isnan(basiccorrect));
    end
end

%% Plotting Superordinate & Basic Category performance (Target trials only)

for i_sub = 1:numSubjects

    % Super - all trials
    super_hitrate(1,i_sub) = nansum(supercorrect_all{1}{i_sub}) / length(supercorrect_all{1}{i_sub});
    super_hitrate(2,i_sub) = nansum(supercorrect_all{2}{i_sub}) / length(supercorrect_all{2}{i_sub});
    super_hitrate(3,i_sub) = nansum(supercorrect_all{3}{i_sub}) / length(supercorrect_all{3}{i_sub});
    % Basic - all trials
    basic_hitrate(1,i_sub) = nansum(basiccorrect_all{1}{i_sub}) / length(basiccorrect_all{1}{i_sub});
    basic_hitrate(2,i_sub) = nansum(basiccorrect_all{2}{i_sub}) / length(basiccorrect_all{2}{i_sub});
    basic_hitrate(3,i_sub) = nansum(basiccorrect_all{3}{i_sub}) / length(basiccorrect_all{3}{i_sub});
    % Basic - supercorrect trials
    basic_supercorr_hitrate(1,i_sub) = nansum(basiccorrect_supercorr{1}{i_sub}) / length(basiccorrect_supercorr{1}{i_sub});
    basic_supercorr_hitrate(2,i_sub) = nansum(basiccorrect_supercorr{2}{i_sub}) / length(basiccorrect_supercorr{2}{i_sub});
    basic_supercorr_hitrate(3,i_sub) = nansum(basiccorrect_supercorr{3}{i_sub}) / length(basiccorrect_supercorr{3}{i_sub});
    
    %Super - undetect
    super_undetectdetect_hitrate(1,i_sub,1) = nansum(supercorrect_undetect{1}{i_sub}) / length(supercorrect_undetect{1}{i_sub});
    super_undetectdetect_hitrate(2,i_sub,1) = nansum(supercorrect_undetect{2}{i_sub}) / length(supercorrect_undetect{2}{i_sub});
    super_undetectdetect_hitrate(3,i_sub,1) = nansum(supercorrect_undetect{3}{i_sub}) / length(supercorrect_undetect{3}{i_sub});
    %Super - detect
    super_undetectdetect_hitrate(1,i_sub,2) = nansum(supercorrect_detect{1}{i_sub}) / length(supercorrect_detect{1}{i_sub});
    super_undetectdetect_hitrate(2,i_sub,2) = nansum(supercorrect_detect{2}{i_sub}) / length(supercorrect_detect{2}{i_sub});
    super_undetectdetect_hitrate(3,i_sub,2) = nansum(supercorrect_detect{3}{i_sub}) / length(supercorrect_detect{3}{i_sub});
    %Super - unrec
    super_unrecrec_hitrate(1,i_sub,1) = nansum(supercorrect_unrec{1}{i_sub}) / length(supercorrect_unrec{1}{i_sub});
    super_unrecrec_hitrate(2,i_sub,1) = nansum(supercorrect_unrec{2}{i_sub}) / length(supercorrect_unrec{2}{i_sub});
    super_unrecrec_hitrate(3,i_sub,1) = nansum(supercorrect_unrec{3}{i_sub}) / length(supercorrect_unrec{3}{i_sub});
    %Super - rec
    super_unrecrec_hitrate(1,i_sub,2) = nansum(supercorrect_rec{1}{i_sub}) / length(supercorrect_rec{1}{i_sub});
    super_unrecrec_hitrate(2,i_sub,2) = nansum(supercorrect_rec{2}{i_sub}) / length(supercorrect_rec{2}{i_sub});
    super_unrecrec_hitrate(3,i_sub,2) = nansum(supercorrect_rec{3}{i_sub}) / length(supercorrect_rec{3}{i_sub});

    %Basic - undetect
    basic_undetectdetect_hitrate(1,i_sub,1) = nansum(basiccorrect_undetect{1}{i_sub}) / length(basiccorrect_undetect{1}{i_sub});
    basic_undetectdetect_hitrate(2,i_sub,1) = nansum(basiccorrect_undetect{2}{i_sub}) / length(basiccorrect_undetect{2}{i_sub});
    basic_undetectdetect_hitrate(3,i_sub,1) = nansum(basiccorrect_undetect{3}{i_sub}) / length(basiccorrect_undetect{3}{i_sub});
    %Basic - detect
    basic_undetectdetect_hitrate(1,i_sub,2) = nansum(basiccorrect_detect{1}{i_sub}) / length(basiccorrect_detect{1}{i_sub});
    basic_undetectdetect_hitrate(2,i_sub,2) = nansum(basiccorrect_detect{2}{i_sub}) / length(basiccorrect_detect{2}{i_sub});
    basic_undetectdetect_hitrate(3,i_sub,2) = nansum(basiccorrect_detect{3}{i_sub}) / length(basiccorrect_detect{3}{i_sub});
    %Basic - unrec
    basic_unrecrec_hitrate(1,i_sub,1) = nansum(basiccorrect_unrec{1}{i_sub}) / length(basiccorrect_unrec{1}{i_sub});
    basic_unrecrec_hitrate(2,i_sub,1) = nansum(basiccorrect_unrec{2}{i_sub}) / length(basiccorrect_unrec{2}{i_sub});
    basic_unrecrec_hitrate(3,i_sub,1) = nansum(basiccorrect_unrec{3}{i_sub}) / length(basiccorrect_unrec{3}{i_sub});
    %Basic - rec
    basic_unrecrec_hitrate(1,i_sub,2) = nansum(basiccorrect_rec{1}{i_sub}) / length(basiccorrect_rec{1}{i_sub});
    basic_unrecrec_hitrate(2,i_sub,2) = nansum(basiccorrect_rec{2}{i_sub}) / length(basiccorrect_rec{2}{i_sub});
    basic_unrecrec_hitrate(3,i_sub,2) = nansum(basiccorrect_rec{3}{i_sub}) / length(basiccorrect_rec{3}{i_sub});
end


titlelist = {'Superordinate','Basic'};
data = {'super_hitrate','basic_hitrate','basic_supercorr_hitrate'};

figure; set(gcf,'Color','w');
hold on

for d = 1:length(data)

disp(data{d});
clearvars individDataPts xval

hitrate_var = eval(data{d});
if strcmp(data{d},'super_hitrate')
    yaxislabel = '% Superordinate Correct';
elseif strcmp(data{d},'basic_hitrate')
    yaxislabel = '% Basic Correct';
elseif strcmp(data{d},'basic_supercorr_hitrate')
    yaxislabel = '% Basic Correct (correct superordinate)';
end

subplot(1,3,d)
hold on

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
      
for r = 1:3               
    xPos = repmat(offsets(r),1,numSubjects);
    individDataPts(r,:) = hitrate_var(r, :)'.*100;
    s = swarmchart(xPos,individDataPts(r,:),30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', condColors{r}); %,'*','LineWidth',2,'Color','k');
    s.XJitter = 'rand';
    s.XJitterWidth = 0.1;
    s.MarkerFaceAlpha = .7;
    xval(r,:) = s.XData;
end
ax=gca;ax.LineWidth=1.2;ax.FontSize = 14; 
ylim([0 100]);

for i_sub = 1:numSubjects
    plot([xval(1,i_sub),xval(2,i_sub),xval(3,i_sub)],[individDataPts(1,i_sub), individDataPts(2,i_sub), individDataPts(3,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
end

for g = 1:3                             
    xPos = repmat(offsets(g),1,numSubjects);
    boxchart(xPos, hitrate_var(g,:).*100,'BoxWidth',0.2,'BoxFaceColor',condColors{g},'JitterOutliers','off','MarkerStyle','none');  
end

set(gca, 'YGrid', 'on');
ylabel(yaxislabel);
yline(50, 'k:','LineWidth',1.5);
legend({cond1name, cond2name, cond3name},'Location','Southeast');
xticklabels({cond1name,cond2name,cond3name});
end

% Testing hit rates against chance level:
% data = {'super_hitrate','basic_hitrate','basic_supercorr_hitrate'};
plotlist{1} = super_hitrate';
plotlist{2} = basic_hitrate';
plotlist{3} = basic_supercorr_hitrate';

for p = 1:length(plotlist)
    temp = plotlist{p};
    for i = 1:3
        data = temp(:,i);
        [pval(p,i), h(p,i), ~] = signrank(data, 50);
    end
end
pval_corr = mafdr(pval(:),'BHFDR','true');
disp(pval_corr);
%% Splitting trials into undetect/detect, or unrec/rec (not included in MS)

data = {'super_undetectdetect_hitrate','basic_undetectdetect_hitrate','super_unrecrec_hitrate', 'basic_unrecrec_hitrate'};

for d = 1:length(data)

disp(data{d});
clearvars individDataPts xval

hitrate_var = eval(data{d});
if contains(data{d},'super')
    yaxislabel = '% Superordinate Correct';
elseif contains(data{d},'basic')
    yaxislabel = '% Basic Correct';
end

if contains(data{d},'detect')
    xtickvar = [{''},{''},{'Undetect'},{''},{'Detect'},{''}];
elseif contains(data{d},'rec')
    xtickvar = [{''},{''},{'Unrec'},{''},{'Rec'},{''}];
end


figure; set(gcf,'Color','w');
hold on

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
groupLabels = {'PAS 0','PAS 1','PAS 2','PAS 3'};

for g = 1:2                  
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        individDataPts{g}(r,:) = hitrate_var(r, :, g)'.*100;
        s = swarmchart(xPos,individDataPts{g}(r,:),30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', condColors{r}); %,'*','LineWidth',2,'Color','k');
        s.XJitter = 'rand';
        s.XJitterWidth = 0.1;
        s.MarkerFaceAlpha = .7;
        xval{g}(r,:) = s.XData;
    end
end
ax=gca;ax.LineWidth=1.2;ax.FontSize = 14; 
ylim([0 100]);

for g = 1:2
    for i_sub = 1:numSubjects
        plot([xval{g}(1,i_sub),xval{g}(2,i_sub),xval{g}(3,i_sub)],[individDataPts{g}(1,i_sub), individDataPts{g}(2,i_sub), individDataPts{g}(3,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
    end
end

for g = 1:2                 
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        boxchart(xPos, hitrate_var(r, :, g).*100,'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off');  
        pval(g,r) = signrank(hitrate_var(r, :, g),50);
    end
end
disp(pval);

set(gca, 'YGrid', 'on');
ylabel(yaxislabel);
yline(50, 'k:','LineWidth',1.5);
legend({cond1name, cond2name, cond3name},'Location','Southeast');
xticklabels(xtickvar);

end
