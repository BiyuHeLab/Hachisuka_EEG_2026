%% Master Script for Behavior analysis (Experiment 3)
% Master script for analysis of behavioral data (EEG behavior)
% Plots Super and Basic Category performance for each condition.
% Plots PAS ratings for each condition.

%Last updated June 20 2025, Ayaka Hachisuka (ahachisu@gmail.com)

clear
close all

clc
load('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/Behavioral/BehavData_30Subs.mat');
numSubjects = 30;
cond1 = 'Normal';
cond2 = 'Masked';
cond3 = 'LSF';
cond_list = {cond1,cond2,cond3};

AllorSuperCorr = NaN; %if NaN, Super Corr trials only. If 0, all trials

if AllorSuperCorr == 0
    filetype = 'AllTrials';
else
    filetype = 'SuperCorrTrialsonly';
end

condColors = {[0.9290, 0.6940, 0.1250],[0.6350, 0.0780, 0.1840],[0.4660 0.6740 0.1880]};

%% 1) Super/Basic % correct by PAS
levelist = {'Basic','Super'};

for ii = 1:2
    level = levelist{ii};
for i_sub = 1:numSubjects
    counter = 0;
    pas_ssub{i_sub} = [];
    pas_ssub_supercorr{i_sub} = [];
    condition_ssub{i_sub} = [];
    condition_ssub_supercorr{i_sub} = [];
    for b = 1:length(BehavData{i_sub}.BlockData)
        pasresp = BehavData{i_sub}.BlockData{b}.ResponseData.PASRespKey;
        pas_ssub{i_sub} = [pas_ssub{i_sub}; pasresp];
        pas_supercorr = BehavData{i_sub}.BlockData{b}.ResponseData.PASRespKey(BehavData{i_sub}.BlockData{b}.ResponseData.SuperRespCorrectNum == 1);
        pas_ssub_supercorr{i_sub} = [pas_ssub_supercorr{i_sub}; pas_supercorr];
        condition = BehavData{i_sub}.BlockData{b}.TrialData.Condition;
        condition_ssub{i_sub} = [condition_ssub{i_sub}; condition];
        condition_supercorr = BehavData{i_sub}.BlockData{b}.TrialData.Condition(BehavData{i_sub}.BlockData{b}.ResponseData.SuperRespCorrectNum == 1);
        condition_ssub_supercorr{i_sub} = [condition_ssub_supercorr{i_sub}; condition_supercorr];

        if strcmp(level,'Super')
            supercorrect = BehavData{i_sub}.BlockData{b}.ResponseData.SuperRespCorrectNum;
        elseif strcmp(level,'Basic')
            supercorrect = BehavData{i_sub}.BlockData{b}.ResponseData.BasicRespCorrectNum;
            supercorrect(supercorrect<0) = AllorSuperCorr;
        end

        pas0corrN = nan(length(condition),1); pas1corrN = nan(length(condition),1);  
        pas2corrN = nan(length(condition),1); pas3corrN = nan(length(condition),1);
        pas0corrM = nan(length(condition),1); pas1corrM = nan(length(condition),1);  
        pas2corrM = nan(length(condition),1); pas3corrM = nan(length(condition),1);
        pas0corrL = nan(length(condition),1); pas1corrL = nan(length(condition),1);  
        pas2corrL = nan(length(condition),1); pas3corrL = nan(length(condition),1);
        for t = 1:length(condition)
            if strcmp(condition(t),'Normal')
                if pasresp(t) == 0 || isnan(pasresp(t))
                    pas0corrN(t) = supercorrect(t);
                elseif pasresp(t) == 1
                    pas1corrN(t) = supercorrect(t);
                elseif pasresp(t) == 2
                    pas2corrN(t) = supercorrect(t);
                elseif pasresp(t) == 3
                    pas3corrN(t) = supercorrect(t);
                end
            elseif strcmp(condition(t),'Masked')
                if pasresp(t) == 0 || isnan(pasresp(t))
                    pas0corrM(t) = supercorrect(t);
                elseif pasresp(t) == 1
                    pas1corrM(t) = supercorrect(t);
                elseif pasresp(t) == 2
                    pas2corrM(t) = supercorrect(t);
                elseif pasresp(t) == 3
                    pas3corrM(t) = supercorrect(t);
                end
            elseif strcmp(condition(t),'LSF')
                if pasresp(t) == 0 || isnan(pasresp(t))
                    pas0corrL(t) = supercorrect(t);
                elseif pasresp(t) == 1
                    pas1corrL(t) = supercorrect(t);
                elseif pasresp(t) == 2
                    pas2corrL(t) = supercorrect(t);
                elseif pasresp(t) == 3
                    pas3corrL(t) = supercorrect(t);
                end
            end
        end
    end

    %remove negatives, which are for incorrect supers -- shouldn't be a
    %problem here, but just in case.
    pas0corrN(pas0corrN<0) = AllorSuperCorr;
    pas1corrN(pas1corrN<0) = AllorSuperCorr;
    pas2corrN(pas2corrN<0) = AllorSuperCorr;
    pas0corrN(pas3corrN<0) = AllorSuperCorr;

    pas0corrM(pas0corrM<0) = AllorSuperCorr;
    pas1corrM(pas1corrM<0) = AllorSuperCorr;
    pas2corrM(pas2corrM<0) = AllorSuperCorr;
    pas0corrM(pas3corrM<0) = AllorSuperCorr;

    pas0corrL(pas0corrL<0) = AllorSuperCorr;
    pas1corrL(pas1corrL<0) = AllorSuperCorr;
    pas2corrL(pas2corrL<0) = AllorSuperCorr;
    pas0corrL(pas3corrL<0) = AllorSuperCorr;

    pas0corrN2(i_sub,:) = nanmean(pas0corrN);
    pas1corrN2(i_sub,:) = nanmean(pas1corrN);
    pas2corrN2(i_sub,:) = nanmean(pas2corrN);
    pas3corrN2(i_sub,:) = nanmean(pas3corrN);

    pas0corrM2(i_sub,:) = nanmean(pas0corrM);
    pas1corrM2(i_sub,:) = nanmean(pas1corrM);
    pas2corrM2(i_sub,:) = nanmean(pas2corrM);
    pas3corrM2(i_sub,:) = nanmean(pas3corrM);

    pas0corrL2(i_sub,:) = nanmean(pas0corrL);
    pas1corrL2(i_sub,:) = nanmean(pas1corrL);
    pas2corrL2(i_sub,:) = nanmean(pas2corrL);
    pas3corrL2(i_sub,:) = nanmean(pas3corrL);

    corrbyPAS_bytrial{ii}(i_sub,:,1,1) = pas0corrN;
    corrbyPAS_bytrial{ii}(i_sub,:,1,2) = pas1corrN;
    corrbyPAS_bytrial{ii}(i_sub,:,1,3) = pas2corrN;
    corrbyPAS_bytrial{ii}(i_sub,:,1,4) = pas3corrN;
    
    corrbyPAS_bytrial{ii}(i_sub,:,2,1) = pas0corrM;
    corrbyPAS_bytrial{ii}(i_sub,:,2,2) = pas1corrM;
    corrbyPAS_bytrial{ii}(i_sub,:,2,3) = pas2corrM;
    corrbyPAS_bytrial{ii}(i_sub,:,2,4) = pas3corrM;
    
    corrbyPAS_bytrial{ii}(i_sub,:,3,1) = pas0corrL;
    corrbyPAS_bytrial{ii}(i_sub,:,3,2) = pas1corrL;
    corrbyPAS_bytrial{ii}(i_sub,:,3,3) = pas2corrL;
    corrbyPAS_bytrial{ii}(i_sub,:,3,4) = pas3corrL;

end

corrbyPAS{ii}(:,1,1) = pas0corrN2;
corrbyPAS{ii}(:,1,2) = pas1corrN2;
corrbyPAS{ii}(:,1,3) = pas2corrN2;
corrbyPAS{ii}(:,1,4) = pas3corrN2;

corrbyPAS{ii}(:,2,1) = pas0corrM2;
corrbyPAS{ii}(:,2,2) = pas1corrM2;
corrbyPAS{ii}(:,2,3) = pas2corrM2;
corrbyPAS{ii}(:,2,4) = pas3corrM2;

corrbyPAS{ii}(:,3,1) = pas0corrL2;
corrbyPAS{ii}(:,3,2) = pas1corrL2;
corrbyPAS{ii}(:,3,3) = pas2corrL2;
corrbyPAS{ii}(:,3,4) = pas3corrL2;

end
%% plot % Correct (Super/Basic) by PAS

ylabelname{1} = '% Super Correct'; %Superordinate d prime';
ylabelname{2} = '% Basic Correct'; %'Basic d prime';

for superorbasic = 1:2

figure; set(gcf,'Color','w');
hold on

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
groupLabels = {'PAS 0','PAS 1','PAS 2','PAS 3'};

for g = 1:4                  
    for r = 1:3              
        xPos = repmat(g + offsets(r),1,numSubjects);
        boxchart(xPos, corrbyPAS{superorbasic}(:,r, g).*100,'BoxWidth',0.2,'BoxFaceColor',condColors{r}, ...
            'JitterOutliers','off','MarkerStyle','none');  

        xPos = repmat(g + offsets(r),1,numSubjects);
        individDataPts = corrbyPAS{superorbasic}(:,r, g)'.*100;
        s = swarmchart(xPos,individDataPts,30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', condColors{r}); %,'*','LineWidth',2,'Color','k');
        s.XJitter = 'rand';
        s.XJitterWidth = 0.1;
        s.MarkerFaceAlpha = .7;
    end
end

xticklabels([{''},{'0'},{'1'},{'2'},{'3'}]);
ylabel(ylabelname{superorbasic});
xlabel('PAS');
set(gca, 'YGrid', 'on');
yline(50,'--');
set(gca,'FontSize',14);
hold off
end

%% Super or Basic Correct by condition

titlelist = {'Superordinate','Basic'};

for i_sub = 1:numSubjects
    counter = 0;
    pas_ssub{i_sub} = [];
    supercorrvar_ssub{i_sub} = [];
    basiccorrvar_ssub{i_sub} = [];
    super_type_ssub{i_sub} = [];
    basic_type_ssub{i_sub} = [];
    condition_ssub{i_sub} = [];
    pasresp_ssub{i_sub} = [];

    for b = 1:length(BehavData{i_sub}.BlockData)
        pasresp = BehavData{i_sub}.BlockData{b}.ResponseData.PASRespKey;
        condition = BehavData{i_sub}.BlockData{b}.TrialData.Condition;
        super_type = BehavData{i_sub}.BlockData{b}.TrialData.TargetImage.Domain;
        basic_type = BehavData{i_sub}.BlockData{b}.TrialData.TargetImage.Class;

        super_correct_var = BehavData{i_sub}.BlockData{b}.ResponseData.SuperRespCorrectNum;
        basic_correct_var = BehavData{i_sub}.BlockData{b}.ResponseData.BasicRespCorrectNum;

        pasresp = BehavData{i_sub}.BlockData{b}.ResponseData.PASRespKey;

        condition_ssub{i_sub} = [condition_ssub{i_sub}; condition];
        pas_ssub{i_sub} = [pas_ssub{i_sub}; pasresp];
        supercorrvar_ssub{i_sub} = [supercorrvar_ssub{i_sub}; super_correct_var];
        basiccorrvar_ssub{i_sub} = [basiccorrvar_ssub{i_sub}; basic_correct_var];
        pasresp_ssub{i_sub} = [pasresp_ssub{i_sub}; pasresp];

        super_type_ssub{i_sub} = [super_type_ssub{i_sub}; super_type];
        basic_type_ssub{i_sub} = [basic_type_ssub{i_sub}; basic_type];
    end
    
    conditionslist = condition_ssub{i_sub};
    supertypelist = super_type_ssub{i_sub};
    basictypelist = basic_type_ssub{i_sub};
    supercorr = supercorrvar_ssub{i_sub};
    basiccorr = basiccorrvar_ssub{i_sub};
    pas = pasresp_ssub{i_sub};

     for c = 1:3
        cond = cond_list{c};

        hit = sum(double(supercorr==1 & strcmp(conditionslist,cond)));
        falsealarm = sum(double(supercorr==0 & strcmp(conditionslist,cond)));
        dprimevar_super(i_sub,c) = norminv(hit/(hit+falsealarm)) - norminv(falsealarm/(hit+falsealarm));
        hitrate_super(i_sub,c) = hit/(hit+falsealarm);

        hit = sum(double(supercorr==1 & strcmp(conditionslist,cond) & strcmp(supertypelist,'Animate')));
        falsealarm = sum(double(supercorr==0 & strcmp(conditionslist,cond) & strcmp(supertypelist,'Animate')));
        hitrate_super2{1}(i_sub,c) = hit/(hit+falsealarm);

        hit = sum(double(supercorr==1 & strcmp(conditionslist,cond) & strcmp(supertypelist,'Inanimate')));
        falsealarm = sum(double(supercorr==0 & strcmp(conditionslist,cond) & strcmp(supertypelist,'Inanimate')));
        hitrate_super2{2}(i_sub,c) = hit/(hit+falsealarm);

        pas_super{1}(i_sub,c) = nanmean(pas(strcmp(conditionslist,cond) & strcmp(supertypelist,'Animate')));
        pas_super{2}(i_sub,c) = nanmean(pas(strcmp(conditionslist,cond) & strcmp(supertypelist,'Inanimate')));

        pas_basic{1}(i_sub,c) = nanmean(pas(strcmp(conditionslist,cond) & strcmp(basictypelist,'Cat')));
        pas_basic{2}(i_sub,c) = nanmean(pas(strcmp(conditionslist,cond) & strcmp(basictypelist,'Dog')));
        pas_basic{3}(i_sub,c) = nanmean(pas(strcmp(conditionslist,cond) & strcmp(basictypelist,'Car')));
        pas_basic{4}(i_sub,c) = nanmean(pas(strcmp(conditionslist,cond) & strcmp(basictypelist,'Truck')));


        clearvars hit falsealarm
        hit = sum(double(basiccorr==1 & strcmp(conditionslist,cond)));
        falsealarm = sum(double(basiccorr==0 & strcmp(conditionslist,cond)));
        falsealarm2 = sum(double(basiccorr<=0 & strcmp(conditionslist,cond)));

        dprimevar_basic(i_sub,c) = norminv(hit/(hit+falsealarm)) - norminv(falsealarm/(hit+falsealarm));
        hitrate_basic(i_sub,c) = hit/(hit+falsealarm);

        hit = sum(double(basiccorr==1 & strcmp(conditionslist,cond) & strcmp(basictypelist,'Cat')));
        falsealarm = sum(double(basiccorr==0 & strcmp(conditionslist,cond) & strcmp(basictypelist,'Cat')));
        hitrate_basic2{1}(i_sub,c) = hit/(hit+falsealarm);

        hit = sum(double(basiccorr==1 & strcmp(conditionslist,cond) & strcmp(basictypelist,'Dog')));
        falsealarm = sum(double(basiccorr==0 & strcmp(conditionslist,cond) & strcmp(basictypelist,'Dog')));
        hitrate_basic2{2}(i_sub,c) = hit/(hit+falsealarm);

        hit = sum(double(basiccorr==1 & strcmp(conditionslist,cond) & strcmp(basictypelist,'Car')));
        falsealarm = sum(double(basiccorr==0 & strcmp(conditionslist,cond) & strcmp(basictypelist,'Car')));
        hitrate_basic2{3}(i_sub,c) = hit/(hit+falsealarm);

        hit = sum(double(basiccorr==1 & strcmp(conditionslist,cond) & strcmp(basictypelist,'Truck')));
        falsealarm = sum(double(basiccorr==0 & strcmp(conditionslist,cond) & strcmp(basictypelist,'Truck')));
        hitrate_basic2{4}(i_sub,c) = hit/(hit+falsealarm);
     end
end

% plotlist = {dprimevar_super,dprimevar_basic};
plotlist = {hitrate_super.*100,hitrate_basic.*100};

for p = 1:length(plotlist)
    temp = plotlist{p};
    for i = 1:3
        data = temp(:,i);
        [pval(p,i), h(p,i), ~] = signrank(data, 50);
    end
end
pval_corr = mafdr(pval(:),'BHFDR','true');


figure; set(gcf,'Color','w');
hold on

for p = 1:length(plotlist)  
subplot(1,length(plotlist),p); hold on
yline(50, 'k:','LineWidth',1.5);
plotdata = plotlist{p};
offsets = linspace(-0.1, 0.1, 3);  % small horizontal shifts for the 3 rows

% ylabel('Performance (d prime)');
ylabel('Performance (% Correct)');

for c = 1:3               
    xPos = repmat(c + offsets(c),1,numSubjects);
    individDataPts(c,:) = squeeze(plotdata(:, c))';
    s = swarmchart(xPos,individDataPts(c,:),30,'black','o','MarkerEdgeColor','none','MarkerFaceColor', condColors{c}); %,'*','LineWidth',2,'Color','k');
    s.XJitter = 'rand';
    s.XJitterWidth = 0.1;
    s.MarkerFaceAlpha = .7;
    xval(c,:) = s.XData;
end

for c = 1:3                 
    xPos = repmat(c + offsets(c),1,numSubjects);
    boxchart(xPos, squeeze(plotdata(:, c)),'BoxWidth',0.4,'BoxFaceColor',condColors{c},'JitterOutliers','off');  
end

for i_sub = 1:30
    plot([xval(1,i_sub),xval(2,i_sub),xval(3,i_sub)],[individDataPts(1,i_sub), individDataPts(2,i_sub), individDataPts(3,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
end

for c = 1:3                 
    xPos = repmat(c + offsets(c),1,numSubjects);
    boxchart(xPos, squeeze(plotdata(:, c)),'BoxWidth',0.4,'BoxFaceColor',condColors{c},'JitterOutliers','off');  
end

ax=gca;ax.LineWidth=1.2;ax.FontSize = 14; 
% ylim([-0.5 5]);
ylim([10 100]);
title(titlelist{p});
set(gca, 'YGrid', 'on');
% % xticklabels ({'OSM','OSOM-stable','OSOM-moved',''});
xticks([1:5]);
xticklabels([{'Normal'},{'Masked'},{'LSF'},{''}]);
set(gca, 'FontName', 'Arial', 'FontSize', 14)
% legend boxoff  
end
disp(pval_corr);

%% hitrate - by animate category & basic category
plottype = 1; % 1 for Hit rate, 2 for PAS
plotsuperorbasic = 2; %1 for SUPERORDINATE, 2 for BASIC

if plottype == 1 % HIT RATE
    yaxisrange = [10 100];
    ylabelval = 'Performance (% Correct)';
    super_plotlist = [hitrate_super2{1}(:,1).*100, hitrate_super2{2}(:,1).*100, ...
        hitrate_super2{1}(:,2).*100, hitrate_super2{2}(:,2).*100, ...
        hitrate_super2{1}(:,3).*100, hitrate_super2{2}(:,3).*100];
    basic_plotlist = [hitrate_basic2{1}(:,1).*100, hitrate_basic2{2}(:,1).*100, hitrate_basic2{3}(:,1).*100, hitrate_basic2{4}(:,1).*100,...
        hitrate_basic2{1}(:,2).*100, hitrate_basic2{2}(:,2).*100, hitrate_basic2{3}(:,2).*100, hitrate_basic2{4}(:,2).*100,...
        hitrate_basic2{1}(:,3).*100, hitrate_basic2{2}(:,3).*100, hitrate_basic2{3}(:,3).*100, hitrate_basic2{4}(:,3).*100];
elseif plottype == 2 % PAS
    yaxisrange = [0 3];
    ylabelval = 'PAS';
    super_plotlist = [pas_super{1}(:,1), pas_super{2}(:,1), ...
        pas_super{1}(:,2), pas_super{2}(:,2), ...
        pas_super{1}(:,3), pas_super{2}(:,3)];
    basic_plotlist = [pas_basic{1}(:,1), pas_basic{2}(:,1), pas_basic{3}(:,1), pas_basic{4}(:,1),...
        pas_basic{1}(:,2), pas_basic{2}(:,2), pas_basic{3}(:,2), pas_basic{4}(:,2),...
        pas_basic{1}(:,3), pas_basic{2}(:,3), pas_basic{3}(:,3), pas_basic{4}(:,3)];
end

plotlist{1} = super_plotlist;
plotlist{2} = basic_plotlist;

for p = 1:length(plotlist)
    temp = plotlist{p};
    for i = 1:3
        data = temp(:,i);
        [pval(p,i), h(p,i), ~] = signrank(data, 50);
    end
end
pval_corr = mafdr(pval(:),'BHFDR','true');


figure; set(gcf,'Color','w');
hold on

for p = plotsuperorbasic%:length(plotlist)  
% subplot(1,length(plotlist),p); hold on

if p == 1
    condColors = {[0.9290, 0.6940, 0.1250],[0.9290, 0.6940, 0.1250],...
    [0.6350, 0.0780, 0.1840],[0.6350, 0.0780, 0.1840], ...
    [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880]};
elseif p == 2
    condColors = {[0.9290, 0.6940, 0.1250],[0.9290, 0.6940, 0.1250], [0.9290, 0.6940, 0.1250], [0.9290, 0.6940, 0.1250],...
    [0.6350, 0.0780, 0.1840],[0.6350, 0.0780, 0.1840], [0.6350, 0.0780, 0.1840], [0.6350, 0.0780, 0.1840],...
    [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880], [0.4660 0.6740 0.1880]};
end

yline(50, 'k:','LineWidth',1.5);
plotdata = plotlist{p};
offsets = linspace(-0.1, 0.1, size(plotlist{p},2));  % small horizontal shifts for the 3 rows

ylabel(ylabelval);

for c = 1:size(plotlist{p},2)               
    xPos = repmat(c + offsets(c),1,numSubjects);
    individDataPts(c,:) = squeeze(plotdata(:, c))';
    s = swarmchart(xPos,individDataPts(c,:),30,'black','o','MarkerEdgeColor','none','MarkerFaceColor', condColors{c}); %,'*','LineWidth',2,'Color','k');
    s.XJitter = 'rand';
    s.XJitterWidth = 0.1;
    s.MarkerFaceAlpha = .7;
    xval(c,:) = s.XData;
end

for c = 1:size(plotlist{p},2)                    
    xPos = repmat(c + offsets(c),1,numSubjects);
    boxchart(xPos, squeeze(plotdata(:, c)),'BoxWidth',0.4,'BoxFaceColor',condColors{c},'JitterOutliers','off');  
end

if p == 1
    for i_sub = 1:30
        plot([xval(1,i_sub),xval(2,i_sub)],[individDataPts(1,i_sub), individDataPts(2,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(3,i_sub),xval(4,i_sub)],[individDataPts(3,i_sub), individDataPts(4,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(5,i_sub),xval(6,i_sub)],[individDataPts(5,i_sub), individDataPts(6,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
    end
elseif p == 2
    for i_sub = 1:30
        plot([xval(1,i_sub),xval(2,i_sub)],[individDataPts(1,i_sub), individDataPts(2,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(2,i_sub),xval(3,i_sub)],[individDataPts(2,i_sub), individDataPts(3,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(3,i_sub),xval(4,i_sub)],[individDataPts(3,i_sub), individDataPts(4,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(5,i_sub),xval(6,i_sub)],[individDataPts(5,i_sub), individDataPts(6,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(6,i_sub),xval(7,i_sub)],[individDataPts(6,i_sub), individDataPts(7,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(7,i_sub),xval(8,i_sub)],[individDataPts(7,i_sub), individDataPts(8,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(9,i_sub),xval(10,i_sub)],[individDataPts(9,i_sub), individDataPts(10,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(10,i_sub),xval(11,i_sub)],[individDataPts(10,i_sub), individDataPts(11,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
        plot([xval(11,i_sub),xval(12,i_sub)],[individDataPts(11,i_sub), individDataPts(12,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);

    end
end
for c = 1:size(plotlist{p},2)                     
    xPos = repmat(c + offsets(c),1,numSubjects);
    boxchart(xPos, squeeze(plotdata(:, c)),'BoxWidth',0.4,'BoxFaceColor',condColors{c},'JitterOutliers','off');  
end

ax=gca;ax.LineWidth=1.2;ax.FontSize = 14; 
% ylim([-0.5 5]);
ylim(yaxisrange);
title(titlelist{p});
set(gca, 'YGrid', 'on');

if p == 1
    xticks([1:6]);
    xticklabels([{'Animal'},{'Vehicle'},{'Animal'},{'Vehicle'},{'Animal'},{'Vehicle'}]);
elseif p == 2
    xticks([1:12]);
    xticklabels([{'Cat'},{'Dog'},{'Car'},{'Truck'},{'Cat'},{'Dog'},{'Car'},{'Truck'},{'Cat'},{'Dog'},{'Car'},{'Truck'}]);
end

set(gca, 'FontName', 'Arial', 'FontSize', 14) 
end

% corrbyPAS stats:
pval_corr = mafdr(pval_pasbycorr{1}(:),'BHFDR','true');
disp(pval_corr);
pval_corr = mafdr(pval_pasbycorr{2}(:),'BHFDR','true');
disp(pval_corr);

%% Plotting number of trials by PAS, by condition
for i_sub = 1:numSubjects
    pas = pas_ssub_supercorr{i_sub};
    cond = condition_ssub_supercorr{i_sub};
    vioPts(1,i_sub,1) = sum(pas == 0 & strcmp(cond,'Normal'))/length(pas);
    vioPts(2,i_sub,1) = sum(pas == 0 & strcmp(cond,'Masked'))/length(pas);
    vioPts(3,i_sub,1) = sum(pas == 0 & strcmp(cond,'LSF'))/length(pas);

    vioPts(1,i_sub,2) = sum(pas == 1 & strcmp(cond,'Normal'))/length(pas);
    vioPts(2,i_sub,2) = sum(pas == 1 & strcmp(cond,'Masked'))/length(pas);
    vioPts(3,i_sub,2) = sum(pas == 1 & strcmp(cond,'LSF'))/length(pas);

    vioPts(1,i_sub,3) = sum(pas == 2 & strcmp(cond,'Normal'))/length(pas);
    vioPts(2,i_sub,3) = sum(pas == 2 & strcmp(cond,'Masked'))/length(pas);
    vioPts(3,i_sub,3) = sum(pas == 2 & strcmp(cond,'LSF'))/length(pas);

    vioPts(1,i_sub,4) = sum(pas == 3 & strcmp(cond,'Normal'))/length(pas);
    vioPts(2,i_sub,4) = sum(pas == 3 & strcmp(cond,'Masked'))/length(pas);
    vioPts(3,i_sub,4) = sum(pas == 3 & strcmp(cond,'LSF'))/length(pas);
end
vioPts = vioPts.*100;

figure; set(gcf,'Color','w');
hold on

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
groupLabels = {'PAS 0','PAS 1','PAS 2','PAS 3'};

for g = 1:4                  
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        boxchart(xPos, vioPts(r, :, g),'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off','MarkerStyle','none');  
    end
end

xticklabels([{''},{'0'},{'1'},{'2'},{'3'}]);
ylabel('% Trials');
xlabel('PAS');

for g = 1:4                   
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        individDataPts = vioPts(r, :, g)';
        s = swarmchart(xPos,individDataPts,30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', condColors{r}); %,'*','LineWidth',2,'Color','k');
        s.XJitter = 'rand';
        s.XJitterWidth = 0.1;
        s.MarkerFaceAlpha = .7;
    end
end
ax=gca;ax.LineWidth=1.2;
% ylim([0 150]);
% title(trialtype);
set(gca, 'YGrid', 'on');
legend({'Normal','Masked','LSF'})
set(gca, 'FontName', 'Arial', 'FontSize', 14)
legend boxoff  

%% Plotting number of trials by PAS, ALL conditions
for i_sub = 1:numSubjects
    pas = pas_ssub_supercorr{i_sub};
    vioPts(i_sub,1) = sum(pas == 0)/length(pas);
    vioPts(i_sub,2) = sum(pas == 1)/length(pas);
    vioPts(i_sub,3) = sum(pas == 2)/length(pas);
    vioPts(i_sub,4) = sum(pas == 3)/length(pas);
end
vioPts = vioPts.*100;

figure; set(gcf,'Color','w');
hold on

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
groupLabels = {'PAS 0','PAS 1','PAS 2','PAS 3'};

for g = 1:4                               
        xPos = repmat(g,1,numSubjects);
        boxchart(xPos, vioPts(:, g),'BoxWidth',0.2,'BoxFaceColor','k','JitterOutliers','off','MarkerStyle','none');  
end

xticklabels([{''},{'0'},{'1'},{'2'},{'3'}]);
ylabel('% Trials');
xlabel('PAS');

for g = 1:4                                
    xPos = repmat(g,1,numSubjects);
    individDataPts = vioPts(:, g)';
    s = swarmchart(xPos,individDataPts,30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', 'k'); %,'*','LineWidth',2,'Color','k');
    s.XJitter = 'rand';
    s.XJitterWidth = 0.1;
    s.MarkerFaceAlpha = .7;
end
ax=gca;ax.LineWidth=1.2;
% ylim([0 150]);
% title(trialtype);
set(gca, 'YGrid', 'on');
legend({'Normal','Masked','LSF'})
set(gca, 'FontName', 'Arial', 'FontSize', 14)
legend boxoff  
%% Percent Trials for Detect/Undetect:
clearvars vioPts

for i_sub = 1:numSubjects
    pas = pas_ssub{i_sub};
    cond = condition_ssub{i_sub};
    vioPts(1,i_sub,1) = sum(pas == 0 & strcmp(cond,'Normal')) / length(pas);
    vioPts(2,i_sub,1) = sum(pas == 0 & strcmp(cond,'Masked')) / length(pas);
    vioPts(3,i_sub,1) = sum(pas == 0 & strcmp(cond,'LSF')) / length(pas);

    vioPts(1,i_sub,2) = sum(pas == 1 | pas == 2 | pas == 3 & strcmp(cond,'Normal')) / length(pas);
    vioPts(2,i_sub,2) = sum(pas == 1 | pas == 2 | pas == 3 & strcmp(cond,'Masked')) / length(pas);
    vioPts(3,i_sub,2) = sum(pas == 1 | pas == 2 | pas == 3 & strcmp(cond,'LSF')) / length(pas);
end

figure; set(gcf,'Color','w');
hold on

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
groupLabels = {'Undetect','Detect'};

for g = 1:2                  
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        boxchart(xPos, vioPts(r, :, g).*100,'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off');  
    end
end

xticklabels([{''},{''},{'Undetect'},{''},{'Detect'},{''}]);
ylabel('% Trials');

for g = 1:2                 
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        individDataPts = vioPts(r, :, g)';
        s = swarmchart(xPos,individDataPts.*100,30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', condColors{r}); %,'*','LineWidth',2,'Color','k');
        s.XJitter = 'rand';
        s.XJitterWidth = 0.1;
        s.MarkerFaceAlpha = .7;
    end
end
ax=gca;ax.LineWidth=1.2;ax.FontSize = 14; 
% ylim([0 150]);
% title(trialtype);
set(gca, 'YGrid', 'on');
legend({'Normal','Masked','LSF'},'Location','Northwest');

%% Percent Trials for Rec/Unrec:
clearvars vioPts

for i_sub = 1:numSubjects
    pas = pas_ssub{i_sub};
    cond = condition_ssub{i_sub};
    vioPts(1,i_sub,1) = sum(pas == 0 | pas == 1 & strcmp(cond,'Normal')) / length(pas);
    vioPts(2,i_sub,1) = sum(pas == 0 | pas == 1 & strcmp(cond,'Masked')) / length(pas);
    vioPts(3,i_sub,1) = sum(pas == 0 | pas == 1 & strcmp(cond,'LSF')) / length(pas);

    vioPts(1,i_sub,2) = sum(pas == 2 | pas == 3 & strcmp(cond,'Normal')) / length(pas);
    vioPts(2,i_sub,2) = sum(pas == 2 | pas == 3 & strcmp(cond,'Masked')) / length(pas);
    vioPts(3,i_sub,2) = sum(pas == 2 | pas == 3 & strcmp(cond,'LSF')) / length(pas);
end

figure; set(gcf,'Color','w');
hold on

offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
groupLabels = {'Unrec','Rec'};

for g = 1:2                  
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        boxchart(xPos, vioPts(r, :, g).*100,'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off');  
    end
end

xticklabels([{''},{''},{'Unrec'},{''},{'Rec'},{''}]);
ylabel('% Trials');

for g = 1:2                 
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        individDataPts = vioPts(r, :, g)';
        s = swarmchart(xPos,individDataPts.*100,30,'black','o', 'MarkerEdgeColor','none','MarkerFaceColor', condColors{r}); %,'*','LineWidth',2,'Color','k');
        s.XJitter = 'rand';
        s.XJitterWidth = 0.1;
        s.MarkerFaceAlpha = .7;
    end
end
ax=gca;ax.LineWidth=1.2;ax.FontSize = 14; 
set(gca, 'YGrid', 'on');
legend({'Normal','Masked','LSF'},'Location','Northeast');

%% Percent correct by undetect/detect

for superorbasic = 1:2

clearvars allcor individDataPts xval

allcorr(1,:,1) = pas0corrN2.*100;
allcorr(1,:,2) = nanmean([pas1corrN2 pas2corrN2 pas3corrN2],2).*100;

allcorr(2,:,1)  = pas0corrM2.*100;
allcorr(2,:,2)  = nanmean([pas1corrM2 pas2corrM2 pas3corrM2],2).*100;

allcorr(3,:,1)  = pas0corrL2.*100;
allcorr(3,:,2)  = nanmean([pas1corrL2 pas2corrL2 pas3corrL2],2).*100;

figure; set(gcf,'Color','w');
hold on
offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows
groupLabels = {'PAS 0','PAS 1','PAS 2','PAS 3'};

for g = 1:2                 
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        boxchart(xPos, allcorr(r, :, g),'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off','MarkerStyle','none');
        pval(g,r) = signrank(allcorr(r, :, g),50);
    end
end

pval_corr = mafdr(pval(:),'BHFDR','true');
disp(pval_corr);

xticklabels([{''},{''},{'Undetect'},{''},{'Detect'}]);
ylabel(ylabelname{superorbasic});
set(gca, 'YGrid', 'on');

for g = 1:2                  
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        individDataPts{g}(r,:) = allcorr(r, :, g)';
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
        boxchart(xPos, allcorr(r, :, g),'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off','MarkerStyle','none');
    end
end
yline(50, 'k:','LineWidth',2.5);
legend({'Normal','Masked','LSF'},'Location','Southeast');
set(gca, 'YGrid', 'on');
set(gca, 'FontName', 'Arial', 'FontSize', 14)
legend boxoff  

%% Percent correct by unrec/rec

clearvars allcor

allcorr(1,:,1) = nanmean([pas0corrN2 pas1corrN2],2).*100;
allcorr(1,:,2) = nanmean([pas2corrN2 pas3corrN2],2).*100;

allcorr(2,:,1)  = nanmean([pas0corrM2 pas1corrM2],2).*100;
allcorr(2,:,2)  = nanmean([pas2corrM2 pas3corrM2],2).*100;

allcorr(3,:,1)  = nanmean([pas0corrL2 pas1corrL2],2).*100;
allcorr(3,:,2)  = nanmean([pas2corrL2 pas3corrL2],2).*100;

figure; set(gcf,'Color','w');
hold on
offsets = linspace(-0.3, 0.3, 3);  % small horizontal shifts for the 3 rows

for g = 1:2                 
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        boxchart(xPos, allcorr(r, :, g),'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off','MarkerStyle','none');
        pval(g,r) = signrank(allcorr(r, :, g),50);
    end
end
disp(pval);

xticklabels([{''},{''},{'Unrec'},{''},{'Rec'}]);
ylabel(ylabelname{superorbasic});
set(gca, 'YGrid', 'on');

for g = 1:2                  
    for r = 1:3               
        xPos = repmat(g + offsets(r),1,numSubjects);
        individDataPts{g}(r,:) = allcorr(r, :, g)';
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
        boxchart(xPos, allcorr(r, :, g),'BoxWidth',0.2,'BoxFaceColor',condColors{r},'JitterOutliers','off','MarkerStyle','none');
    end
end
yline(50, 'k:','LineWidth',2.5);
legend({'Normal','Masked','LSF'},'Location','Southeast');
set(gca, 'YGrid', 'on');
set(gca, 'FontName', 'Arial', 'FontSize', 14)
legend boxoff  

end

%% MEAN PAS PLOTS:
levelist = {'Super','Basic'};
pas_all = cell(3,1); %[];
pas_superINcorr = cell(3,1);
pas_basicINcorr = cell(3,1);
pas_supercorr = cell(3,1); %[];
pas_basiccorr = cell(3,1); %[];
for cond = 1:3
    pas_all{cond} = cell(numSubjects,1);
    pas_superINcorr{cond} = cell(numSubjects,1);
    pas_basicINcorr{cond} = cell(numSubjects,1);
    pas_supercorr{cond} = cell(numSubjects,1);
    pas_basiccorr{cond} = cell(numSubjects,1);
end

for i_sub = 1:numSubjects
    counter = 0;
    pas_ssub{i_sub} = [];
    condition_ssub{i_sub} = [];

    for b = 1:length(BehavData{i_sub}.BlockData)
        clearvars temp_pas temp_pas_supercorr temp_pas_basiccorr
        pasresp = BehavData{i_sub}.BlockData{b}.ResponseData.PASRespKey;
        % pas_ssub{i_sub} = [pas_ssub{i_sub}; pasresp];
        condition = BehavData{i_sub}.BlockData{b}.TrialData.Condition;
        % condition_ssub{i_sub} = [condition_ssub{i_sub}; condition];

        supercorrect = BehavData{i_sub}.BlockData{b}.ResponseData.SuperRespCorrectNum;
        basiccorrect = BehavData{i_sub}.BlockData{b}.ResponseData.BasicRespCorrectNum;
        basiccorrect(basiccorrect<0)=0;
        temp_pas{1} = pasresp(contains(condition,'Normal'));
        temp_pas_superINcorr{1} = pasresp(supercorrect~=1 & contains(condition,'Normal'));
        temp_pas_supercorr{1} = pasresp(supercorrect==1 & contains(condition,'Normal'));
        temp_pas_basicINcorr{1} = pasresp(basiccorrect~=1 & contains(condition,'Normal'));
        temp_pas_basiccorr{1} = pasresp(basiccorrect==1 & contains(condition,'Normal'));

        temp_pas{2} = pasresp(contains(condition,'Masked'));
        temp_pas_superINcorr{2} = pasresp(supercorrect~=1 & contains(condition,'Masked'));
        temp_pas_supercorr{2} = pasresp(supercorrect==1 & contains(condition,'Masked'));
        temp_pas_basicINcorr{2} = pasresp(basiccorrect~=1 & contains(condition,'Masked'));
        temp_pas_basiccorr{2} = pasresp(basiccorrect==1 & contains(condition,'Masked'));

        temp_pas{3} = pasresp(contains(condition,'LSF'));
        temp_pas_superINcorr{3} = pasresp(supercorrect~=1 & contains(condition,'LSF'));
        temp_pas_supercorr{3} = pasresp(supercorrect==1 & contains(condition,'LSF'));
        temp_pas_basicINcorr{3} = pasresp(basiccorrect~=1 & contains(condition,'LSF'));
        temp_pas_basiccorr{3} = pasresp(basiccorrect==1 & contains(condition,'LSF'));

        pas_all{1}{i_sub} = [pas_all{1}{i_sub}; temp_pas{1}];
        pas_superINcorr{1}{i_sub} = [pas_superINcorr{1}{i_sub}; temp_pas_superINcorr{1}];
        pas_supercorr{1}{i_sub} = [pas_supercorr{1}{i_sub}; temp_pas_supercorr{1}];
        pas_basicINcorr{1}{i_sub} = [pas_basicINcorr{1}{i_sub}; temp_pas_basicINcorr{1}];
        pas_basiccorr{1}{i_sub} = [pas_supercorr{1}{i_sub}; temp_pas_basiccorr{1}];

        pas_all{2}{i_sub} = [pas_all{2}{i_sub}; temp_pas{2}];
        pas_superINcorr{2}{i_sub} = [pas_superINcorr{2}{i_sub}; temp_pas_superINcorr{2}];
        pas_supercorr{2}{i_sub} = [pas_supercorr{2}{i_sub}; temp_pas_supercorr{2}];
        pas_basicINcorr{2}{i_sub} = [pas_basicINcorr{2}{i_sub}; temp_pas_basicINcorr{2}];
        pas_basiccorr{2}{i_sub} = [pas_basiccorr{2}{i_sub}; temp_pas_basiccorr{2}];

        pas_all{3}{i_sub} = [pas_all{3}{i_sub}; temp_pas{3}];
        pas_superINcorr{3}{i_sub} = [pas_superINcorr{3}{i_sub}; temp_pas_superINcorr{3}];
        pas_supercorr{3}{i_sub} = [pas_supercorr{3}{i_sub}; temp_pas_supercorr{3}];
        pas_basicINcorr{3}{i_sub} = [pas_basicINcorr{3}{i_sub}; temp_pas_basicINcorr{3}];
        pas_basiccorr{3}{i_sub} = [pas_basiccorr{3}{i_sub}; temp_pas_basiccorr{3}];
    end
end

for cond = 1:3
    for i_sub = 1:numSubjects
        pas_allsub(i_sub,cond) = nanmean(pas_all{cond}{i_sub});
        pas_allsub_superINcorr(i_sub,cond) = nanmean(pas_superINcorr{cond}{i_sub});
        pas_allsub_supercorr(i_sub,cond) = nanmean(pas_supercorr{cond}{i_sub});
        pas_allsub_basicINcorr(i_sub,cond) = nanmean(pas_basicINcorr{cond}{i_sub});
        pas_allsub_basiccorr(i_sub,cond) = nanmean(pas_basiccorr{cond}{i_sub});
    end
end

%% Plotting PAS ratings by trial condition:

analysistypelist = {'PAS with all trials','PAS with correct super trials','PAS with correct basic trials', ...
    'PAS with INcorrect super trials','PAS with INcorrect basic trials'};
figure; set(gcf,'Color','w');
hold on

for aa = 1:length(analysistypelist)
    analysistype = analysistypelist{aa};

    subplot(2,3,aa); hold on
if contains(analysistype,'PAS with all trials')
    normal_pas = pas_allsub(:,1);
    masked_pas = pas_allsub(:,2);
    lsf_pas = pas_allsub(:,3);
    % save('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Behavior_stats/pas_allsub_alltrials.mat','pas_allsub');
elseif contains(analysistype,'PAS with correct super trials')
    normal_pas = pas_allsub_supercorr(:,1);
    masked_pas = pas_allsub_supercorr(:,2);
    lsf_pas = pas_allsub_supercorr(:,3);
    % save('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Behavior_stats/pas_allsub_supercorrtrials.mat','pas_allsub_supercorr');
elseif contains(analysistype,'PAS with correct basic trials')
    normal_pas = pas_allsub_basiccorr(:,1);
    masked_pas = pas_allsub_basiccorr(:,2);
    lsf_pas = pas_allsub_basiccorr(:,3);
    % save('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Behavior_stats/pas_allsub_basiccorrtrials.mat','pas_allsub_basiccorr');
elseif contains(analysistype,'PAS with INcorrect super trials')
    normal_pas = pas_allsub_superINcorr(:,1);
    masked_pas = pas_allsub_superINcorr(:,2);
    lsf_pas = pas_allsub_superINcorr(:,3);
    % save('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Behavior_stats/pas_allsub_supercorrtrials.mat','pas_allsub_supercorr');
elseif contains(analysistype,'PAS with INcorrect basic trials')
    normal_pas = pas_allsub_basicINcorr(:,1);
    masked_pas = pas_allsub_basicINcorr(:,2);
    lsf_pas = pas_allsub_basicINcorr(:,3);
    % save('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Behavior_stats/pas_allsub_supercorrtrials.mat','pas_allsub_supercorr');
end

all_pas = [normal_pas masked_pas lsf_pas];

condColor = {[0.9290, 0.6940, 0.1250],[0.6350, 0.0780, 0.1840],[0.4660 0.6740 0.1880]};

clearvars individDataPts xval

plotdata = all_pas;
offsets = linspace(-0.1, 0.1, 3);  % small horizontal shifts for the 3 rows

for c = 1:3               
    xPos = repmat(c + offsets(c),1,numSubjects);
    individDataPts(c,:) = squeeze(plotdata(:, c));
    s = swarmchart(xPos,individDataPts(c,:),30,'black','o','MarkerEdgeColor','none','MarkerFaceColor', condColors{c}); %,'*','LineWidth',2,'Color','k');
    s.XJitter = 'rand';
    s.XJitterWidth = 0.1;
    s.MarkerFaceAlpha = .7;
    xval(c,:) = s.XData;
end

for c = 1:3                 
    xPos = repmat(c + offsets(c),1,numSubjects);
    boxchart(xPos, squeeze(plotdata(:, c)),'BoxWidth',0.4,'BoxFaceColor',condColors{c},'JitterOutliers','off');  
end

for i_sub = 1:30
    plot([xval(1,i_sub),xval(2,i_sub),xval(3,i_sub)],[individDataPts(1,i_sub), individDataPts(2,i_sub), individDataPts(3,i_sub)],'-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
end

for c = 1:3                 
    xPos = repmat(c + offsets(c),1,numSubjects);
    boxchart(xPos, squeeze(plotdata(:, c)),'BoxWidth',0.4,'BoxFaceColor',condColors{c},'JitterOutliers','off');  
end

ax=gca;ax.LineWidth=1.2;ax.FontSize = 14; 
title(analysistype);
set(gca, 'YGrid', 'on');
ylim([0 3]);
xticks([1:5]);
xticklabels([{'Normal'},{'Masked'},{'LSF'},{''}]);
set(gca, 'FontName', 'Arial', 'FontSize', 14)
end

