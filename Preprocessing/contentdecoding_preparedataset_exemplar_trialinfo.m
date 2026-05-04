%% This is a general ERP analysis script.
%Last edited, 2/17/23, Ayaka Hachisuka
paths_C2F_EEG = C2F_EEG_SetPaths;
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
vars = vars_C2F_EEG;

supercorrectTrials = 0;
basiccorrectTrials = 1;
recTrials = 0;
unrecTrials = 0;

timewindow = [-0.4 0.8];

%% Add paths:
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Scripts/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Scripts/SubFunctions/');

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts/Preprocessing/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts/Preprocessing/SubFunctions/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts/GeneralFunctions/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/external/mne/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/fileio/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/utilities/');
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fieldtrip-release/plotting');

paths_C2F_EEG = C2F_EEG_SetPaths;
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
vars = vars_C2F_EEG;

%% Load data
% subjList = {'01','02','03','04','07','08','10','11','13','14','15','17','18','20','22','24','25','26','27','28'};
subjList = {'01','02','03','04','07','08','10','11','13','14','15','17','18','20','22','24','25','26','27','28','19','21','23','29','30','31','33','34','35','36','37','38','39'};
% subjList = {'31','33','34','35','36','37','38'};
for s = 1:length(subjList)
    subnum = subjList{s};
    % datall = load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/Clean/FinalClean_original_deletebadchanandtrials/sub',subnum,'_dataEEG_cleanedFinal']);
    datall = load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub',subnum,'_dataEEG_cleanedFinal2']);

numSubjects = length(subjList);

%% separating trials by super or basic categories:
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');

    data = datall.dataEEG_cleanedFinal2;
    subjNumStr = SubLabel_Add0BeforeSingleDigit(s); % add a 0 before single digit subjNums
    % %check time indices for all chanels (to get a sense of trial length:
    % for ch = 1:length(timetemp)
    %     trial_length(ch) = size(data.trial{ch},2);
    % end
    % bar(trial_length);

    %Find trials with ONLY correct super (and exluding 1-1-1 no resp trials)
    validtrial = ones(size(data.trialinfo,1),1); %all trials are valid unless otherwise noted
    numtrial = size(data.trialinfo,1);
    if supercorrectTrials == 1
        for t = 1:numtrial %number of trials
            info = data.trialinfo(t,:);
            if info(15) == 0 || isnan(info(15)) %super is incorrect
                validtrial(t) = 0;
            end
        end
    elseif basiccorrectTrials == 1
        for t = 1:numtrial %number of trials
            info = data.trialinfo(t,:);
            if info(17) == 0 || isnan(info(17)) %super is incorrect
                validtrial(t) = 0;
            end
        end
    elseif recTrials == 1
        for t = 1:numtrial %number of trials
            info = data.trialinfo(t,:);
            if info(19) == 0 || isnan(info(19)) || info(19) == 1 %if UNrecognized, make invalid
                validtrial(t) = 0;
            end
        end
     elseif unrecTrials == 1
        for t = 1:numtrial %number of trials
            info = data.trialinfo(t,:);
            if info(19) == 2 || info(19) == 3 %if recognized, make invalid
                validtrial(t) = 0;
            end
        end
    end  
    validtrialIND = find(validtrial==1);
    
    % ANTI-ALIAS FILTERING AT LP 35 Hz, ON CONTINUOUS DATA:
    cfg = [];
    cfg.channel = {'all', '-Mast_L', '-Mast_R', '-HEOG_L','-HEOG_R', '-VEOG_Inf','-VEOG_Sup',  '-VEOG_Sp', '-ECG_Inf','-ECG_Sup'};  %EEG & EOG only
    cfg.continuous  = 'yes'; %preprocess the entire recording.
    cfg.lpfilter    = 'yes'; %yes' or 'no' (default = 'no')
    cfg.lpfreq      = 35; %scalar value for low pass frequency (there is no default, so needs to be always specified)
    cfg.lpfilttype  = 'firws'; %string, filter type (default is set in ft_preproc_lowpassfilter)
    data = ft_preprocessing(cfg, data);

    cfg = [];
    cfg.trials = validtrialIND;
    
    % EPOCH DATA FROM -0.4 to 0.8s.
    cfg.channel = 'all';
    cfg.latency = timewindow;
    
    data2 = ft_selectdata(cfg,data);
    data2.info.chs = data2.label;
    data2.trialinfo = data.trialinfo(validtrialIND,:);
    
    % RESAMPLE DATA (& BASELINE CORRECT HERE.)
    cfg = [];
    cfg.resamplefs = 200;
    cfg.demean = 'yes';
    cfg.baselinewindow = [-0.4 0];
    data2 = ft_resampledata(cfg, data2);
      
    trialinfo_temp = data2.trialinfo;
    newinfo = zeros(size(trialinfo_temp,1),1);
    for t = 1:size(trialinfo_temp,1)
        if trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 1 %normal, dog1
            newinfo(t) = 1;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 1 %masked, dog1
            newinfo(t) = 2;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 1 %lsf, dog1
            newinfo(t) = 3;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 2 %normal, dog2
            newinfo(t) = 4;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 2 %masked, dog2
            newinfo(t) = 5;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 2 %lsf, dog2
            newinfo(t) = 6;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 3 %normal, dog3
            newinfo(t) = 7;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 3 %masked, dog3
            newinfo(t) = 8;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 3 %lsf, dog3
            newinfo(t) = 9;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 4 %normal, cat1
            newinfo(t) = 10;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 4 %masked, cat1
            newinfo(t) = 11;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 4 %lsf, cat1
            newinfo(t) = 12;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 5 %normal, cat2
            newinfo(t) = 13;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 5 %masked, cat2
            newinfo(t) = 14;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 5 %lsf, cat2
            newinfo(t) = 15;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 6 %normal, cat3
            newinfo(t) = 16;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 6 %masked, cat3
            newinfo(t) = 17;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 6 %lsf, cat3
            newinfo(t) = 18;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 7 %normal, car1
            newinfo(t) = 19;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 7 %masked, car1
            newinfo(t) = 20;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 7 %lsf, car1
            newinfo(t) = 21;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 8 %normal, car2
            newinfo(t) = 22;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 8 %masked, car2
            newinfo(t) = 23;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 8 %lsf, car2
            newinfo(t) = 24;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 9 %normal, car3
            newinfo(t) = 25;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 9 %masked, car3
            newinfo(t) = 26;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 9 %lsf, car3
            newinfo(t) = 27;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 10 %normal, truck1
            newinfo(t) = 28;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 10 %masked, truck1
            newinfo(t) = 29;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 10 %lsf, truck1
            newinfo(t) = 30;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 11 %normal, truck2
            newinfo(t) = 31;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 11 %masked, truck2
            newinfo(t) = 32;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 11 %lsf, truck2
            newinfo(t) = 33;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,12) == 12 %normal, truck3
            newinfo(t) = 34;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,12) == 12 %masked, truck3
            newinfo(t) = 35;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,12) == 12 %lsf, truck3
            newinfo(t) = 36;
        end
    end
    data2.origtrialinfo = data2.trialinfo;
    data2.trialinfo = [];
    data2.trialinfo = newinfo;
    
    cfg = [];
    cfg.keeptrials = 'yes';
    cfg.removemean = 'no';
    data3 = ft_timelockanalysis(cfg,data2);
    data3.origtrialinfo = data2.origtrialinfo;
    if supercorrectTrials == 1 %if you want to do this, do the same for basic correct. 
        save(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_inputv2/','sub',num2str(subjNumStr),'_dataEEG_SuperCorrect'],'data3');
    elseif basiccorrectTrials == 1 %if you want to do this, do the same for basic correct.
        save(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_inputv2/','sub',num2str(subjNumStr),'_dataEEG_BasicCorrect'],'data3');
    elseif recTrials == 1
        save(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_inputv2/','sub',num2str(subjNumStr),'_dataEEG_Rec'],'data3');
    elseif unrecTrials == 1
        save(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_inputv2/','sub',num2str(subjNumStr),'_dataEEG_Unrec'],'data3');
    else
        save(['/isilon/LFMI/VMdrive/Ayaka/EEG/ADAM_data/ContentDecoding/ExemplarDecoding200Hz/final_decoding_inputv2/','sub',num2str(subjNumStr),'_dataEEG_ALL'],'data3');
    end
end

%% sanity check average ERP plot:
% paths_C2F_EEG = C2F_EEG_SetPaths;
% vars = C2F_EEG_SetVars(paths_C2F_EEG);
% vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
% initDir = pwd;
% 
% cfg = [];
% cfg.demean = 'yes';
% cfg.baselinewindow = [-0.5 0];
% data_temp2 = ft_preprocessing(cfg, data2);

cfg = [];
ft_singleplotER(cfg,data3);

