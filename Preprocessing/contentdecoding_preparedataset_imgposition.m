%% This is a general ERP analysis script.
%Last edited, 2/17/23, Ayaka Hachisuka
% paths_C2F_EEG = C2F_EEG_SetPaths;
% vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
% vars = vars_C2F_EEG;

correctTrials = 1;
recTrials = 0;
unrecTrials = 0;
undetectTrials = 0;

timewindow = [-0.4 0.8];

%% Load data

% if ~exist('datall', 'var')
%     datall{1} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub01_dataEEG_cleanedFinal2');
%     datall{2} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub02_dataEEG_cleanedFinal2');
%     datall{3} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub03_dataEEG_cleanedFinal2');
%     datall{4} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub04_dataEEG_cleanedFinal2');
%     datall{5} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub07_dataEEG_cleanedFinal2');
%     datall{6} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub08_dataEEG_cleanedFinal2');
%     datall{7} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub10_dataEEG_cleanedFinal2');
%     datall{8} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub11_dataEEG_cleanedFinal2');
%     datall{9} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub13_dataEEG_cleanedFinal2');
%     datall{10} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub14_dataEEG_cleanedFinal2');
%     datall{11} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub15_dataEEG_cleanedFinal2');
%     datall{12} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub17_dataEEG_cleanedFinal2');
%     datall{13} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub18_dataEEG_cleanedFinal2');
%     datall{14} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub20_dataEEG_cleanedFinal2');
%     datall{15} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub22_dataEEG_cleanedFinal2');
%     datall{16} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub24_dataEEG_cleanedFinal2');
%     datall{17} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub25_dataEEG_cleanedFinal2');
%     datall{18} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub26_dataEEG_cleanedFinal2');
%     datall{19} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub27_dataEEG_cleanedFinal2');
%     datall{20} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub28_dataEEG_cleanedFinal2');
%     datall{21} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub19_dataEEG_cleanedFinal2');
%     datall{22} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub23_dataEEG_cleanedFinal2');
%     datall{23} = load('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub29_dataEEG_cleanedFinal2');
% end

addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Scripts/');
addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Scripts/SubFunctions/');

addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/fromThomas/Scripts/Preprocessing/');
addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/fromThomas/Scripts/Preprocessing/SubFunctions/');
addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/fromThomas/Scripts/GeneralFunctions/');
addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/fieldtrip-release/');
addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/fieldtrip-release/external/mne/');
addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/fieldtrip-release/fileio/');
addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/fieldtrip-release/utilities/');
addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/fieldtrip-release/plotting');
paths_C2F_EEG = C2F_EEG_SetPaths;
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
vars = vars_C2F_EEG;
%% separating trials by super or basic categories:
%% Load data
subjList = {'01','02','03','04','07','08','10','11','13','14','15','17','18','20','22','24','25','26','27','28','19','21','23','29','30','31','33','34','35','36','37','38'};

for s = 1:length(subjList)
    subnum = subjList{s};
    data = load(['/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Data/EEG_preproc/Clean/FinalClean_repairedchan2_keepchan/sub',subnum,'_dataEEG_cleanedFinal2']).dataEEG_cleanedFinal2;
    addpath('/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/Preprocessing');
    subjNumStr = SubLabel_Add0BeforeSingleDigit(s); % add a 0 before single digit subjNums
    % %check time indices for all chanels (to get a sense of trial length:
    % for ch = 1:length(timetemp)
    %     trial_length(ch) = size(data.trial{ch},2);
    % end
    % bar(trial_length);

    %Find trials with ONLY correct super (and exluding 1-1-1 no resp trials)
    validtrial = ones(size(data.trialinfo,1),1); %all trials are valid unless otherwise noted
    numtrial = size(data.trialinfo,1);
    if correctTrials == 1
        for t = 1:numtrial %number of trials
            info = data.trialinfo(t,:);
            if isnan(info(15)) ||  info(15) == 0 %super is incorrect
                validtrial(t) = 0;
            end
        end
    elseif recTrials == 1
        for t = 1:numtrial %number of trials
            info = data.trialinfo(t,:);
            if info(19) == 0 || isnan(info(19)) || info(19) == 1  %if UNrecognized, make invalid
                validtrial(t) = 0;
            end
        end
     elseif unrecTrials == 1
        for t = 1:numtrial %number of trials
            info = data.trialinfo(t,:);
            if  info(19) == 2 || info(19) == 3 %if recognized, make invalid
                validtrial(t) = 0;
            end
        end
     elseif undetectTrials == 1
        for t = 1:numtrial %number of trials
            info = data.trialinfo(t,:);
            if  info(19) == 1 || info(19) == 2 || info(19) == 3 %if recognized, make invalid
                validtrial(t) = 0;
            end
        end
    % else
    %     for t = 1:numtrial %number of trials
    %         info = data.trialinfo(t,:);
    %         if  isnan(info(19)) %if no PAS response, make invalid
    %             validtrial(t) = 0;
    %         end
    %     end
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
        if trialinfo_temp(t,9) == 1 && trialinfo_temp(t,14) == 1 %normal, top
            newinfo(t) = 1;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,14) == 1 %masked, top
            newinfo(t) = 2;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,14) == 1 %lsf, top
            newinfo(t) = 3;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,14) == 2 %normal, bottom
            newinfo(t) = 4;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,14) == 2 %masked, bottom
            newinfo(t) = 5;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,14) == 2 %lsf, bottom
            newinfo(t) = 6;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,14) == 3 %normal, left
            newinfo(t) = 7;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,14) == 3 %masked, left
            newinfo(t) = 8;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,14) == 3 %lsf, left
            newinfo(t) = 9;
        elseif trialinfo_temp(t,9) == 1 && trialinfo_temp(t,14) == 4 %normal, right
            newinfo(t) = 10;
        elseif trialinfo_temp(t,9) == 2 && trialinfo_temp(t,14) == 4 %masked, right
            newinfo(t) = 11;
        elseif trialinfo_temp(t,9) == 3 && trialinfo_temp(t,14) == 4 %lsf, right
            newinfo(t) = 12;
        end
    end
    data2.trialinfo = [];
    data2.trialinfo = newinfo;

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
    data2.trialinfo(:,2) = newinfo;
    data2.origtrialinfo = trialinfo_temp;
    
    cfg = [];
    cfg.keeptrials = 'yes';
    cfg.removemean = 'no';

    data3 = ft_timelockanalysis(cfg,data2);
    data3.origtrialinfo = trialinfo_temp;
    if correctTrials == 1 %if you want to do this, do the same for basic correct. But this will probably remain unused, not bothering with it now.
        save(['/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/ADAM_data/PositionDecoding/final_decoding_input/','sub',num2str(subjNumStr),'_dataEEG_SuperCorrect'],'data3');
    elseif recTrials == 1
        save(['/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/ADAM_data/PositionDecoding/final_decoding_input/','sub',num2str(subjNumStr),'_dataEEG_Rec2'],'data3');
    elseif unrecTrials == 1
        save(['/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/ADAM_data/PositionDecoding/final_decoding_input/','sub',num2str(subjNumStr),'_dataEEG_UNRec2'],'data3');
     elseif undetectTrials == 1
        save(['/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/ADAM_data/PositionDecoding/final_decoding_input/','sub',num2str(subjNumStr),'_dataEEG_UNdetect'],'data3');
    else
        save(['/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/ADAM_data/PositionDecoding/final_decoding_input/','sub',num2str(subjNumStr),'_dataEEG_ALL'],'data3');
        % save(['/isilon/LFMI/VMdrive/Ayaka/C2F_eeg/ADAM_data/PositionDecoding/final_decoding_input/','sub',num2str(subjNumStr),'_dataEEG_ALL_withPASresp'],'data3');
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

