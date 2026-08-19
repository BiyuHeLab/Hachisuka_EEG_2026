%% PLOT group analysis for time-generalization matrices with stats:
clear;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');
paths_EEG = EEG_SetPaths;
vars_EEG = EEG_SetVars(paths_EEG);
vars = vars_EEG;

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
config=cosmo_config();

conditionList = {'Normal','Masked','LSF'};
superorbasicList = {'Super'};%,'Basic'};

n = 256;
condColors = {[0.9290, 0.6940, 0.1250],[0.6350, 0.0780, 0.1840],[0.4660 0.6740 0.1880]};

for c = 1:length(conditionList)
condition = conditionList{c};
for sorb = 1:length(superorbasicList)
        superorbasic = superorbasicList{sorb};
num_subjects = 33;
subjectslist = setdiff(1:num_subjects,[22 32 33]);

%load data
if contains(superorbasic, 'Basic')
    addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts/Preprocessing');
    filedir1 = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/Time_generalization_matrices/Timeradius_is_1/AllTrials/Basic1/', condition,'/']);
    % filedir1 = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/Time_generalization_matrices/Basic1/',condition,'/']);
    % data1 = load([filedir1,'/allsub_dgm_ds_cell'],'dgm_ds_cell').dgm_ds_cell;
    % groupds1 = load([filedir1,'/allsub_group_ds'],'group_ds').group_ds;
    scounter = 0;
    for k = subjectslist
        scounter = scounter + 1;
        subjNumStr = SubLabel_Add0BeforeSingleDigit(k);
        load([filedir1,'sub',num2str(subjNumStr),'_',condition,'_cdt_ds.mat']);
        cdt_ds.samples = cdt_ds.samples';
        [data, labels, values]=cosmo_unflatten(cdt_ds,1);
        dgm_ds_cell1{scounter} = cdt_ds;
        ds=cosmo_dim_transpose(cdt_ds,{'train_time', 'test_time' },2);
        % for one-sample t-test
        ds.sa.targets=1;

        % each participant is independent
        ds.sa.chunks=scounter;
        ds.a.fdim.values{1,1} = round(ds.a.fdim.values{1,1},3); %round time variable, because there is a miniscule jitter in the "0" value, creating problems with data stacking across subs.
        ds.a.fdim.values{2,1} = round(ds.a.fdim.values{2,1},3);
        data1{scounter}=ds;
    end
    groupds1=cosmo_stack(data1);

    % filedir2 = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/Time_generalization_matrices/Basic2/',condition,'/']);
    filedir2 = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/Time_generalization_matrices/Timeradius_is_1/AllTrials/Basic2/', condition,'/']);
    % data2 = load([filedir1,'/allsub_dgm_ds_cell'],'dgm_ds_cell').dgm_ds_cell;
    % groupds2 = load([filedir2,'/allsub_group_ds'],'group_ds').group_ds;
    scounter = 0;
    for k = subjectslist
        scounter = scounter + 1;
        subjNumStr = SubLabel_Add0BeforeSingleDigit(k);
        load([filedir2,'sub',num2str(subjNumStr),'_',condition,'_cdt_ds.mat']);
        cdt_ds.samples = cdt_ds.samples';
        [data, labels, values]=cosmo_unflatten(cdt_ds,1);
        dgm_ds_cell2{scounter} = cdt_ds;
        ds=cosmo_dim_transpose(cdt_ds,{'train_time', 'test_time' },2);
        % for one-sample t-test
        ds.sa.targets=1;

        % each participant is independent
        ds.sa.chunks=scounter;
        ds.a.fdim.values{1,1} = round(ds.a.fdim.values{1,1},3); %round time variable, because there is a miniscule jitter in the "0" value, creating problems with data stacking across subs.
        ds.a.fdim.values{2,1} = round(ds.a.fdim.values{2,1},3);
        data2{scounter}=ds;
    end
    groupds2=cosmo_stack(data2);

    for s = 1:length(subjectslist)
        tempdgm{1,s}.samples = (dgm_ds_cell1{1,s}.samples + dgm_ds_cell2{1,s}.samples) / 2;
        tempdgm{1,s}.a = dgm_ds_cell1{1,s}.a;
        tempdgm{1,s}.sa = dgm_ds_cell1{1,s}.sa;
    end
    dgm_ds_cell = tempdgm;

    samp(:,:,1) = groupds1.samples;
    samp(:,:,2) = groupds2.samples;
    tempds.samples = mean(samp,3);
    tempds.a = groupds1.a;
    tempds.fa = groupds1.fa;
    tempds.sa = groupds1.sa;
    group_ds = tempds;
else
    addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts/Preprocessing');
    filedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/Time_generalization_matrices/Timeradius_is_1/AllTrials/Super/', condition,'/']);
    % filedir = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan200Hz_balexemp/Time_generalization_matrices/',superorbasic,'/',condition,'/']);
    scounter = 0;
    for k = subjectslist
        scounter = scounter + 1;
        subjNumStr = SubLabel_Add0BeforeSingleDigit(k);
        load([filedir,'sub',num2str(subjNumStr),'_',condition,'_cdt_ds.mat']);
        cdt_ds.samples = cdt_ds.samples';
        [data, labels, values]=cosmo_unflatten(cdt_ds,1);
        dgm_ds_cell{scounter} = cdt_ds;
        ds=cosmo_dim_transpose(cdt_ds,{'train_time', 'test_time' },2);
        % for one-sample t-test
        ds.sa.targets=1;

        % each participant is independent
        ds.sa.chunks=scounter;
        ds.a.fdim.values{1,1} = round(ds.a.fdim.values{1,1},3); %round time variable, because there is a miniscule jitter in the "0" value, creating problems with data stacking across subs.
        ds.a.fdim.values{2,1} = round(ds.a.fdim.values{2,1},3);
        group_cell{scounter}=ds;
    end
    group_ds=cosmo_stack(group_cell);

    % load([filedir,'/allsub_dgm_ds_cell'],'dgm_ds_cell');
    % load([filedir,'/allsub_group_ds'],'group_ds');
end

%%
% run multiple comparison correction
nbrhood=cosmo_cluster_neighborhood(group_ds);
ds_result = cosmo_montecarlo_cluster_stat(group_ds,nbrhood,'h0_mean',0.5,'niter',1000,'cluster_stat','maxsize','threshold',1,'p_uncorrected',0.05);
% save(fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan/Time_generalization_matrices/',superorbasic,'/',condition,'/ds_result']),'ds_result');

% addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/examples/redblue/');

% for easier unflattenening, the time dimensions are moved back
% from feature to sample dimensions
ds_result_time = cosmo_dim_transpose(ds_result,...
                                    { 'train_time', 'test_time' });

% flatten into a 2D array
[arr,dim_labels,dim_values]=cosmo_unflatten(ds_result_time,1);
arr(arr>1.64) = 1;
arr(arr~=1) = 0;

%% average across all subs:
subjNum = length(dgm_ds_cell);
for k = 1:subjNum
    cdt_ds = dgm_ds_cell{k};
    % cdt_ds.fa.test_time = cdt_ds.fa.test_time';
    % cdt_ds.fa.train_time = cdt_ds.fa.train_time';
    [data, labels, values]=cosmo_unflatten(cdt_ds,1);
    data_allsubs(:,:,k) = data;
end

% [data, labels, values]=cosmo_unflatten(group_ds,1);
meandata = mean(data_allsubs,3);

n = 256;
if c == 1
    targetcolor = condColors{1};
elseif c == 2
    targetcolor = condColors{2};
elseif c == 3
    targetcolor = condColors{3};
end

figure; set(gcf,'Color','w');
hold on

img = imagesc(meandata);
axis equal

% Outline sig elements
contour(arr, 1, 'LineColor', 'k', 'LineWidth', 1.5);

axis equal tight;
hold off;
ytick=1:20:numel(dim_values{1});
ylabel(strrep(dim_labels{1},'_',' '));
set(gca,'Ytick',ytick,'YTickLabel',round(dim_values{1}(ytick),2));

xtick=1:20:numel(dim_values{2});
xlabel(strrep(dim_labels{2},'_',' '));
set(gca,'Xtick',xtick,'XTickLabel',round(dim_values{2}(xtick),2));
xtickangle(45);

axisval = dim_values{2};
zeroval = find(axisval == 0);
lastval = length(axisval);
line([zeroval, zeroval], [0, lastval], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');
line([0, lastval], [zeroval, zeroval], 'LineStyle', '--', 'LineWidth', 1, 'Color', 'k');

cmap = [linspace(1,targetcolor(1),n)', ...
        linspace(1,targetcolor(2),n)', ...
        linspace(1,targetcolor(3),n)'];
colormap(cmap)

grid on

ax = gca;
ax.Layer = 'top';
ax.GridColor = [0 0 0];
ax.GridAlpha = 0.15;
grid on

colorbar();
title([condition, '-', superorbasic, ' (N = ', num2str(subjNum),')']);

clim([0.495 0.53]);

end
end

% fullpath = fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan/Time_generalization_matrices/',condition,'_',superorbasic,'_timegenALL']);
% saveas(gcf, fullpath);
