%% Position decoding
%normal, top - 1
%masked, top - 2
%lsf, top - 3
%normal, bottom - 4
%masked, bottom - 5
%lsf, bottom - 6
%normal, left - 7
%masked, left - 8
%lsf, left - 9 
%normal, right - 10
%masked, right - 11
%lsf, right - 12

conditionlist = {'Masked','LSF'};
% conditionlist = {'Normal'};
for cc = 1:length(conditionlist)

condition = conditionlist{cc};

run_pca = 0;
time_radius = 0;        
trial_bin_num = 4;
trial_resample = 1;
nchunks = 10;
nminval = 1;
num_subjects = 32;

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing')
paths_EEG = EEG_SetPaths;
vars_EEG = EEG_SetVars(paths_EEG);
vars = vars_EEG;

%% get timelock data in CoSMoMVPA format
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/external/');

% set configuration
config=cosmo_config();
config.data_path = fullfile('/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/PositionDecoding/final_decoding_input/');
% save(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/PositionDecoding/final_decoding_input/','sub',num2str(subjNumStr),'_dataEEG_ALL'],'data3');

data_path=fullfile(config.data_path);
index2label = {'Position1','Position2'};

% reset citation list
cosmo_check_external('-tic');

for s = 1:num_subjects

    subjNum = s;
    if length(num2str(subjNum)) == 1 || ...
            (mod(subjNum, 1) ~= 0 && floor(subjNum) < 10)
        subjNumStr = ['0' num2str(subjNum)];
    else
        subjNumStr = num2str(subjNum);
    end

for posvalue1 = 1:4 %one vs other schematic
    for posvalue2 = posvalue1:4 %one vs other schematic
        if posvalue1~=posvalue2
            disp(sprintf(['Position value: ',num2str(posvalue1),' vs. ', num2str(posvalue2)]));
            data_fn=fullfile([data_path,'sub', subjNumStr, '_dataEEG_ALL.mat']);
        data_tl1=load(data_fn);
        data_tl = data_tl1.data3;
    
        %only take the conditions of interest:
        if contains(condition, 'Normal')
            cond1 = 1; cond2 = 4; cond3 = 7; cond4 = 10;
        elseif contains(condition, 'Masked')
            cond1 = 2; cond2 = 5; cond3 = 8; cond4 = 11;
        elseif contains(condition, 'LSF')
            cond1 = 3; cond2 = 6; cond3 = 9; cond4 = 12;
        end
    
        indices = find(data_tl.trialinfo ~= cond1 & data_tl.trialinfo ~= cond2 & data_tl.trialinfo ~= cond3 & data_tl.trialinfo ~= cond4);
        data_tl.trialinfo(indices) = [];
        data_tl.trial(indices,:,:) = [];
        data_tl.sampleinfo(indices,:) = [];
    
        fourcondvals = unique(data_tl.trialinfo);
        data_tl.trialinfo(data_tl.trialinfo == fourcondvals(1)) = 1;
        data_tl.trialinfo(data_tl.trialinfo == fourcondvals(2)) = 2;
        data_tl.trialinfo(data_tl.trialinfo == fourcondvals(3)) = 3;
        data_tl.trialinfo(data_tl.trialinfo == fourcondvals(4)) = 4;
    
        indices = find(data_tl.trialinfo ~= posvalue1 & data_tl.trialinfo ~= posvalue2);
        data_tl.trialinfo(indices) = [];
        data_tl.trial(indices,:,:) = [];
        data_tl.sampleinfo(indices,:) = [];
    
        data_tl.trialinfo2(data_tl.trialinfo == posvalue1) = 1;
        data_tl.trialinfo2(data_tl.trialinfo == posvalue2) = 2; % set data elements equal to othervalues to 2
    
        data_tl.trialinfo = data_tl.trialinfo2';
        % convert to cosmomvpa struct
        addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/toolboxes/CoSMoMVPA/mvpa/');
        ds_tl=cosmo_meeg_dataset(data_tl);

        %% Prepare MVPA
        % set the target (trial condition)
        ds_tl.sa.targets=data_tl.trialinfo; % 1=pre, 2=post
        
        % in addition give a label to each trial
        ds_tl.sa.labels=cellfun(@(x)index2label(x),num2cell(ds_tl.sa.targets));
        
        % just to check everything is ok
        cosmo_check_dataset(ds_tl);
        
        % reset chunks: use 12 chunks
        ds_tl.sa.chunks=(1:size(ds_tl.samples,1))';
        ds_tl.sa.chunks=cosmo_chunkize(ds_tl,nchunks);
        
        %%%%%%%%
        discardsub = 0;
        target1 = find(ds_tl.sa.targets == 1);
        target2 = find(ds_tl.sa.targets == 2);
        
        %4-fold CV, at least 2 in each target x 4 = 8 minimum. Discard if 9 or
        %less.
        if length(target1) < trial_bin_num || length(target2) < trial_bin_num
            discardsub = 1;
        end
        
        if discardsub > 0
            disp('Not enough trials, discarding subject.')
            continue;
        else
        average_train_args_cell={{'average_train_count',trial_bin_num,...            %  { 2
                                        'average_train_resamplings',trial_resample}}; %  { 2
        n_average_train_args=numel(average_train_args_cell);
        
        chantypes = {'eeg'};
        
        nchantypes=length(chantypes);
        
        partitions=cosmo_nchoosek_partitioner(ds_tl,1);
        partitions=cosmo_balance_partitions(partitions, ds_tl,'nmin',nminval);
        
        npartitions=numel(partitions);
        fprintf('There are %d partitions\n', numel(partitions.train_indices));
        fprintf('# train samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                                partitions.train_indices)));
        fprintf('# test samples:%s\n', sprintf(' %d', cellfun(@numel, ...
                                                partitions.test_indices)));
        
        measure=@cosmo_crossvalidation_measure;
        measure_args=struct();
        measure_args.classifier=@cosmo_classify_lda3;
        measure_args.partitions=partitions;
        measure_args.output = 'balanced_accuracy';
        measure_args.normalize = 'zscore';
        if run_pca == 1
            measure_args.pca_explained_ratio = 0.95;
        end
        if ~isempty(average_train_args_cell)
            average_train_args=average_train_args_cell{1};
            measure_args=cosmo_structjoin(measure_args, average_train_args);
        end
    
        nbrhood=cosmo_interval_neighborhood(ds_tl,'time',...
                                                'radius',time_radius);
        
        sl_map=cosmo_searchlight(ds_tl,nbrhood,measure,measure_args);
        
        for t = 1:length(sl_map.activation_pattern)
            for fold = 1:length(sl_map.activation_pattern{t})
                temp_map(:,:,fold) = sl_map.activation_pattern{t}{fold};
            end
            act_map(:,t) = mean(temp_map,3);
        end

        sl_map.activation_pattern = act_map;

        fprintf('The output has feature dimensions: %s\n', ...
                        cosmo_strjoin(sl_map.a.fdim.labels,', '));
        
        f = figure('visible','off');
        hold on
        time_values=sl_map.a.fdim.values{1}; % first dim (channels got nuked)
        plot(time_values,sl_map.samples,'LineWidth',1.5,'Color','b');
        xline(0, '--', 'LineWidth', 1, 'Color', 'k');
        yline(0.5, '--', 'LineWidth', 1, 'Color', 'k');
        xlim([min(time_values),max(time_values)]);
        ylabel('Balanced Accuracy (chance=.5)');
        
        xlabel('time');
        
        descr=sprintf(['Subject ', num2str(s),', Position value ', num2str(posvalue1), ' vs. ', num2str(posvalue2)]);
        title(descr);
    
        % save fig
        fig_save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan/PositionDecoding1vs1/',condition,'/noPCA_35Hz_100Hz_', num2str(nchunks), 'chunks_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_position',num2str(posvalue1),num2str(posvalue2),'_LDAba_AllTrials/figs/']);
        if ~exist(fig_save_dir, 'dir'); mkdir(fig_save_dir); end
        filename = sprintf('sub%d.jpg',s);
        fullpath = fullfile(fig_save_dir, filename);
        saveas(gcf, fullpath);
    
        % Save data output files:
        save_dir=fullfile(['/isilon/LFMI/VMdrive/Ayaka/EEG/CosmoMVPA_results/AllChan/PositionDecoding1vs1/',condition,'/noPCA_35Hz_100Hz_', num2str(nchunks), 'chunks_', num2str(time_radius), 'timerad_', num2str(trial_bin_num), 'trialresample_morechan_position',num2str(posvalue1),num2str(posvalue2),'_LDAba_AllTrials/']);
    
        if ~exist(save_dir, 'dir'); mkdir(save_dir); end
        save([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_mvpaoutput'],'sl_map','-v7.3');
        %% 
        % save([save_dir,'sub',num2str(subjNumStr),'_dataEEG_',condition,'_ds_tl'],'ds_tl','-v7.3');
        end
    end
end
end
end
end
