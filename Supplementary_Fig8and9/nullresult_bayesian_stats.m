%% Computing Bayesian statistics for C2F: Supplementary Fig. 6 and 7
% Supplementary Fig. 6 is for mutlivariate timeseries decoding analysis
% with undetected trials.
% Supplemetnary Fig. 7 is for CCGP and PS analysis.

% Last updated July 17 2025, Ayaka Hachisuka (ahachisu@gmail.com)

clear;

%% Loading data for CCGP & PS (Supplementary Fig 6)

addpath(genpath('C:\Users\ahach\OneDrive\Desktop\C2F_manuscript\BFF_repo-master\codes\local_functions'));

cd ('C:\Users\ahach\OneDrive\Desktop\C2F_manuscript\BF_stats\CCGP_alltrials');
% cd('C:\Users\ahach\OneDrive\Desktop\C2F_manuscript\BF_stats\BasicControl_alltrials');
load('data_for_bf_ccgp.mat');
load('data_for_bf_control.mat');
load('data_for_ps.mat');

for c = 1:3
    x1 = data_for_bf_within{c} - 0.5;
    bayes_within{c} = bayesfactor_R_wrapper(x1',...
        'returnindex',2,'verbose',1,'args','mu=0,rscale="medium",nullInterval=c(-Inf,0.5)');
    x2 = data_for_bf_across{c} - 0.5;
    bayes_across{c} = bayesfactor_R_wrapper(x2',...
        'returnindex',2,'verbose',1,'args','mu=0,rscale="medium",nullInterval=c(-Inf,0.5)');
    bayes_ps{c} = bayesfactor_R_wrapper(ps{c}',...
        'returnindex',2,'verbose',1,'args','mu=0,rscale="medium",nullInterval = c(-0.1, 0.1)');
end

%% Plotting BF for CCGP
condColor{1} = [155 155 155]./255; %[82 202 255]./255;
condColor{2} = [0 18 154]./255;

conditionList = {'Normal','Masked','LSF'};
x = linspace(-0.4,0.8,241);
figure; set(gcf,'Color','w'); hold on;
for c = 1:3
    subplot(3,1,c); hold on;

    bf1 = bayes_within{c}';
    logBF01_within = log(bf1);

    bf2 = bayes_across{c}';
    logBF01_basic = log(bf2);

    stem(x,logBF01_within,'Marker','o','Color',.6*[1 1 1],'BaseValue',0,'MarkerSize',5,'MarkerFaceColor',condColor{1},'Clipping','off');
    plot(x,logBF01_within,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',condColor{1},'Clipping','off');
    stem(x,logBF01_basic,'Marker','o','Color',.6*[1 1 1],'BaseValue',0,'MarkerSize',5,'MarkerFaceColor',condColor{2},'Clipping','off');
    plot(x,logBF01_basic,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',condColor{2},'Clipping','off');
    xlabel('Time (ms)')
    ylabel('Log BF_{10}')   
    xlim([-0.4 0.8]);
    ylim([-15 10]);
    yline(log10(3), 'k');
    yline(-log10(3), 'k');
    title(conditionList{c});
    grid on
end

%% Plotting BF for PS
condColors = {[0.9290, 0.6940, 0.1250],[0.6350, 0.0780, 0.1840],[0.4660 0.6740 0.1880]};

figure; set(gcf,'Color','w'); hold on;
for c = 1:3
    subplot(3,1,c); hold on;

    bf1 = bayes_ps{c}';
    logbf = log(bf1);

    stem(x,logbf,'Marker','o','Color',.6*[1 1 1],'BaseValue',0,'MarkerSize',5,'MarkerFaceColor','r','Clipping','off');
    plot(x,logbf,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',condColors{c},'Clipping','off');
    xlabel('Time (ms)')
    ylabel('Log BF_{10}')
    xlim([-0.4 0.8]);
    ylim([-4 6]);
    yline(log10(3), 'k');
    yline(-log10(3), 'k');
    title(conditionList{c});
    grid on
end

%% Loading data for undetected Super and Basic (Supplemetnary Fig. 7):

clear;
cd ('C:\Users\ahach\OneDrive\Desktop\C2F_manuscript\BF_stats\');
load('data_undetect_superandbasic.mat');
%load('data_for_bf_clmm.mat');
data_for_bf = data_undetect_superandbasic;

x1 = data_for_bf{1} - 0.5;
bayes{1} = bayesfactor_R_wrapper(x1',...
    'returnindex',2,'verbose',1,'args','mu=0,rscale="medium",nullInterval=c(-Inf,0.5)');
x2 = data_for_bf{2} - 0.5;
bayes{2} = bayesfactor_R_wrapper(x2',...
    'returnindex',2,'verbose',1,'args','mu=0,rscale="medium",nullInterval=c(-Inf,0.5)');

%% Plotting BF for undetected Super and Basic decoding outputs:

condColor{1} = [50 130 246]./255;
condColor{2} = [240 134 80]./255;

conditionList = {'Super','Basic'};
x = linspace(-0.4,0.8,241);
figure; set(gcf,'Color','w'); hold on;

bf1 = bayes{1}';
logBF01_super = log(bf1);

bf2 = bayes{2}';
logBF01_basic = log(bf2);

subplot(2,1,1); hold on;
stem(x,logBF01_super,'Marker','o','Color',.6*[1 1 1],'BaseValue',0,'MarkerSize',5,'MarkerFaceColor',condColor{1},'Clipping','off');
plot(x,logBF01_super,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',condColor{1},'Clipping','off');

xlabel('Time (ms)')
ylabel('Log BF_{10}')   
xlim([-0.4 0.8]);
ylim([-15 10]);
yline(log10(3), 'k');
yline(-log10(3), 'k');
title(conditionList{1});
grid on

subplot(2,1,2); hold on;
stem(x,logBF01_basic,'Marker','o','Color',.6*[1 1 1],'BaseValue',0,'MarkerSize',5,'MarkerFaceColor',condColor{2},'Clipping','off');
plot(x,logBF01_basic,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',condColor{2},'Clipping','off');
xlabel('Time (ms)')
ylabel('Log BF_{10}')   
xlim([-0.4 0.8]);
ylim([-15 10]);
yline(log10(3), 'k');
yline(-log10(3), 'k');
title(conditionList{2});
grid on

