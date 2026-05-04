%% Preprocess EEG data
% This script re-integrates trial information that got accidentally deleted
% during previous preprocessing steps. Rely on absolute trial number
% preserved in trialinfo.
        
%NOTES:
    %6 Missing electrodes (used as EOG/ECG): O9,Ol1h, Ol2h, O10, P9, P10)
    %i.e., 128-6-2(mastoids) = 120 active EEG electrodes in an actiCAP 128 channel Standard2 config
    
%% Preparation:
%Header
clear; clc

global subjNum

addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/') %main project path
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts/') %all scripts
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts/Preprocessing/SubFunctions/') %subfunctions
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts/GeneralFunctions/') %general functions
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');

% C2F_Setup_setVars
% C2F_Setup_setAllDirectories;

paths_C2F_EEG = C2F_EEG_SetPaths;
vars = C2F_EEG_SetVars(paths_C2F_EEG);
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
% cd /isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts

initDir = pwd;

ft_default.trackusage = 'no';
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');
% subjNum = vars.validSubjs; %Specify current subject

for subjNum = 37%[1 2 3 4 7 8 10 11 13 14 15 17 18 19 20 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38]
% subjNum = 1;
subjNumStr = SubLabel_Add0BeforeSingleDigit(subjNum); % add a 0 before single digit subjNums

%1. Load in behavioral data
% load([vars.C2F_behaviorDir, 'C2F_behavData.mat']);
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');
rawBehavData = C2F_EEG_createBehavData(subjNum);

rawBehavData2(subjNum) = rawBehavData;
rawBehavData = rawBehavData2;

% Load in EEG data
taskTypeNum = 2; %because only 2(Experiment) has EEG data
taskType = vars.Task.TaskTypes{taskTypeNum};

%Specify EEG data path
subjDir = ['sub' subjNumStr '/Experiment_data/' taskType  '/'];
rawDataFileName = ['C2F_EEG_00' subjNumStr '.vhdr'] ; % ['sub' subjNumStr 'pos_allBlocks.vhdr'];
preProcessingSettingsName = 'default_EEG';
eegDataDir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/RawData/', subjDir,'EEG_rawData/'];
cd(eegDataDir)
RawEEGDataFileName = [eegDataDir rawDataFileName];    
EEGTriggerFileName = [rawDataFileName(1:end-4) 'vmrk'];

% Preprocessing begins here:

% subjNumStr = '18b';

cfg = [];
cfg.dataset = RawEEGDataFileName;
cfg.continuous  = 'yes'; %Tell FT that the data is one continuous recording
% cfg.trialfun = 'C2F_trialfun_defineblocks';
cfg.trialfun = 'C2F_trialfun_definetrials';
cfg = ft_definetrial(cfg);
cfg.hdr = ft_read_header(cfg.dataset);

%2.1 Read out trial-specific (i.e., condition) and behavioral (i.e., response) data and append it to the trialinfo subfield,
%determine fit by checking/comparing the start times of the respective trial in the EEG and behavioral file
for block = 1:length(rawBehavData(subjNum).Experiment.blockData)
    for trial_perblock = 1:length(rawBehavData(subjNum).Experiment.blockData(block).trialData.trialnumber)
        
        trial_starttime = rawBehavData(subjNum).Experiment.blockData(block).triggerData.StartTime_block...
                                   +rawBehavData(subjNum).Experiment.blockData(block).triggerData.StartTimeAligned2BlockStart_startphase(trial_perblock);
        trial = find(cfg.trl(:,1) == trial_starttime);   
        
        % trial_starttimeall(block,trial_perblock) = trial_starttime;
        if isempty(trial)
            emptytrials{block,trial_perblock} = 1;
        elseif ~isempty(trial)
                %Condition               
                if strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.condition(trial_perblock),'Normal')                      
                cfg.trl(trial,12) = 1;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.condition(trial_perblock),'Masked')   
                cfg.trl(trial,12) = 2;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.condition(trial_perblock),'LSF') 
                cfg.trl(trial,12) = 3;
                end
                
                %Domain
                if strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.domain(trial_perblock),'Animate')                      
                cfg.trl(trial,13) = 1;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.domain(trial_perblock),'Inanimate')   
                cfg.trl(trial,13) = 2;
                end 
                
                %Class                
                if strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.class(trial_perblock),'Dog')                      
                cfg.trl(trial,14) = 1;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.class(trial_perblock),'Cat')   
                cfg.trl(trial,14) = 2;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.class(trial_perblock),'Car') 
                cfg.trl(trial,14) = 3;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.class(trial_perblock),'Truck') 
                cfg.trl(trial,14) = 4;
                end
                
                %Image Exemplar                
                if ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Dog1'))                     
                cfg.trl(trial,15) = 1;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Dog2'))   
                cfg.trl(trial,15) = 2;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Dog3')) 
                cfg.trl(trial,15) = 3;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Cat1'))                      
                cfg.trl(trial,15) = 4;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Cat2'))   
                cfg.trl(trial,15) = 5;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Cat3')) 
                cfg.trl(trial,15) = 6;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Car1'))                      
                cfg.trl(trial,15) = 7;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Car2'))  
                cfg.trl(trial,15) = 8;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Car3')) 
                cfg.trl(trial,15) = 9;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Truck1'))                      
                cfg.trl(trial,15) = 10;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Truck2'))   
                cfg.trl(trial,15) = 11;
                elseif ~cellfun(@isempty,strfind(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage(trial_perblock),'Truck3')) 
                cfg.trl(trial,15) = 12;
                end                
                
                %Contrast Level
                cfg.trl(trial,16) = rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage_contrastlvl(trial_perblock);
                
                %Image Location
                if strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage_location(trial_perblock),'Top')                      
                cfg.trl(trial,17) = 1;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage_location(trial_perblock),'Bottom')   
                cfg.trl(trial,17) = 2;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage_location(trial_perblock),'Left') 
                cfg.trl(trial,17) = 3;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.targetimage_location(trial_perblock),'Right') 
                cfg.trl(trial,17) = 4;
                end
                
                %Response correctness and RT
                cfg.trl(trial,18) = rawBehavData(subjNum).Experiment.blockData(block).responseData.SuperordinateRespCorrectnessNum(trial_perblock);
                cfg.trl(trial,19) = rawBehavData(subjNum).Experiment.blockData(block).responseData.SuperordinateRTinMS(trial_perblock);        
                cfg.trl(trial,20) = rawBehavData(subjNum).Experiment.blockData(block).responseData.BasicRespCorrectnessNum(trial_perblock);
                cfg.trl(trial,21) = rawBehavData(subjNum).Experiment.blockData(block).responseData.BasicRTinMS(trial_perblock);  
                cfg.trl(trial,22) = rawBehavData(subjNum).Experiment.blockData(block).responseData.PASResp(trial_perblock);
                cfg.trl(trial,23) = rawBehavData(subjNum).Experiment.blockData(block).responseData.PASRTinMS(trial_perblock);
                
                %Epoch durations
                cfg.trl(trial,24) = rawBehavData(subjNum).Experiment.blockData(block).trialData.duration_startphase(trial_perblock);
                cfg.trl(trial,25) = rawBehavData(subjNum).Experiment.blockData(block).trialData.duration_prestimjitter(trial_perblock);
                cfg.trl(trial,26) = rawBehavData(subjNum).Experiment.blockData(block).trialData.duration_targetimage(trial_perblock);
                cfg.trl(trial,27) = rawBehavData(subjNum).Experiment.blockData(block).trialData.duration_delaymask(trial_perblock);
                cfg.trl(trial,28) = rawBehavData(subjNum).Experiment.blockData(block).trialData.duration_poststimjitter(trial_perblock);
        end
    end
end

% % If there are misisng blocks, delete them
% for i = 1:size(cfg.trl,1)
%     if nansum(cfg.trl(i,25:28)) == 0
%         cfg.trl(i,:) = nan;
%     else
%         cfg.trl(i,:) = cfg.trl(i,:);
%     end
% end
% 
% cfg.trl = cfg.trl(~all(isnan(cfg.trl), 2), :);

trl_struct = cfg.trl; %store trialstructure in separate variable for later.

%for subject 4, delete trials from the last block -- DEPRECATED, THIS
%PROBLEM IS SOLVED UPSTREAM! There was just a small fluke block with no
%real behavior data.
    % if subjNum == 4
    %     for i = 1:size(trl_struct,1)
    %         if sum(trl_struct(i,12:28)) == 0
    %             trl_struct(i,:) = NaN;
    %         end
    %     end
    %     allnanrows = all(isnan(trl_struct), 2);
    %     trl_struct(allnanrows,:) = [];
    % end

% save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/trl_struct_files/sub',num2str(subjNumStr),'_trl_struct'],'trl_struct', '-v7.3');

%% Load final data & fill with a complete trl_struct matrix.

load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean_repairedchan_keepchan/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal2'],'dataEEG_cleanedFinal2');

%make sure only EEG channels make it to the last dataset:
cfg  = [];
cfg.channel = {'all', '-Mast_L', '-Mast_R', '-HEOG_L','-HEOG_R', '-VEOG_Inf','-VEOG_Sup',  '-VEOG_Sp', '-ECG_Inf','-ECG_Sup'}; 
dataEEG_cleanedFinal2 = ft_preprocessing(cfg, dataEEG_cleanedFinal2);

final_trialnums = dataEEG_cleanedFinal2.trialinfo(:,8);

% if subjNum == 4
%     behavtrials = trl_struct(:,11);
%     eegtrials = final_trialnums;
%     for e = 1:size(final_trialnums,1)
%         if ismember(final_trialnums(e),behavtrials)
%             final_trlstruct(e,12:28) = trl_struct(e,12:28);
%         else
%             dataEEG_cleanedFinal2.sampleinfo(e,:) = NaN;
%             dataEEG_cleanedFinal2.trial{e} = NaN;
%             dataEEG_cleanedFinal2.time{e} = NaN;
%         end
%     end
%     dataEEG_cleanedFinal2.sampleinfo(any(isnan(dataEEG_cleanedFinal2.sampleinfo), 2), :) = [];
%     dataEEG_cleanedFinal2.trial = dataEEG_cleanedFinal2.trial(~cellfun(@(x) isequaln(x, NaN), dataEEG_cleanedFinal2.trial));
%     dataEEG_cleanedFinal2.time = dataEEG_cleanedFinal2.time(~cellfun(@(x) isequaln(x, NaN), dataEEG_cleanedFinal2.time));
if subjNum == 27 || subjNum == 18 %if the trial numbering restarts after two sessions are combined (subj 18, 27)
    temp = find(final_trialnums == 1);
    secondstart = temp(2); %index of the trial number at the beginning of session 2
    firstsessTrialNum = final_trialnums(secondstart-1); %trial number at the end of session 1
    final_trialnums2 = [final_trialnums(1:secondstart-1); final_trialnums(secondstart:end) + firstsessTrialNum];
    final_trialrows = ismember(trl_struct(:,11), final_trialnums2); 
    final_trlstruct = trl_struct(final_trialrows,:);
else
    final_trialrows = ismember(trl_struct(:,11), final_trialnums); 
    final_trlstruct = trl_struct(final_trialrows,:);
end


%delete first 3 columns
final_trlstruct2 = final_trlstruct; %first save original
final_trlstruct(:,1:3) = []; %delete first 3 columns

dataEEG_cleanedFinal2.trialinfo = final_trlstruct;
% dataEEG_cleanedFinal2.trialinfo2 = final_trlstruct2;

save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean_repairedchan2_keepchan/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal2'],'dataEEG_cleanedFinal2', '-v7.3');
end
%% Test again by plotting ERPs

% temporalInfo_Artifacts.manual_postICA = dataEEG_cleanedFinal.cfg.previous.artfctdef;
% temporalInfo_Artifacts.manual_postICA.channel = badchan;
% save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/','sub',num2str(subjNumStr),'_SummaryRejectedArtifacts'],'temporalInfo_Artifacts');

%% Test plot 1: ERP
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal'],'dataEEG_cleanedFinal');
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/trl_struct_files/sub',num2str(subjNumStr),'_trl_struct'],'trl_struct');
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/SemiautoClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedSemiauto'],'dataAll_EEGcleanedSemiauto');

% figure;
% hold on
for s = 37
subjNumStr = SubLabel_Add0BeforeSingleDigit(s); % add a 0 before single digit subjNums

load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean_repairedchan2_keepchan/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal2'],'dataEEG_cleanedFinal2');
load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/trl_struct_files/sub',num2str(subjNumStr),'_trl_struct'],'trl_struct');
datatemp = dataEEG_cleanedFinal2;
goodtrials = datatemp.cfg.trials;
% datatemp = dataAll_EEGcleanedSemiauto;
cfg = [];
% cfg.trl = trl_struct(goodtrials,:);
cfg.trl = trl_struct;
% cfg.channel = {'Oz'};
% cfg.channel = {'all', '-HEOG_L','-HEOG_R', '-VEOG_Inf','-VEOG_Sp','-VEOG_Sup', '-ECG_Inf','-ECG_Sup'}; %select only EEG channels
cfg.channel = {'all'};
cfg.demean = 'yes';
cfg.baselinewindow = [-0.5 0]; %min duration of prestim fixation point relativ to target image
data_temp2 = ft_preprocessing(cfg, datatemp);

cfg.latency = [-0.5000 0.5000];
dataEEG_cleanedFinal2 = ft_selectdata(cfg, data_temp2);

% dataEEG_cleanedFinal2.trialinfo = final_trlstruct;

Normaltrials = find(dataEEG_cleanedFinal2.trialinfo(:,9) == 1);
Maskedtrials = find(dataEEG_cleanedFinal2.trialinfo(:,9) == 2);
LSFtrials = find(dataEEG_cleanedFinal2.trialinfo(:,9) == 3);

cfg = [];
cfg.latency = [-.5 0.5];
cfg.vartrllength = 2;
cfg.trials = Normaltrials;
test_normal = ft_timelockanalysis(cfg,dataEEG_cleanedFinal2);
cfg.trials = Maskedtrials;
test_masked = ft_timelockanalysis(cfg,dataEEG_cleanedFinal2);
cfg.trials = LSFtrials;
test_LSF = ft_timelockanalysis(cfg,dataEEG_cleanedFinal2);

% subplot(5,7,s)
cfg = [];
ft_singleplotER(cfg,test_normal,test_masked,test_LSF);
        legend('Normal','Masked','LSF','Location','northwest')
title(s);
end
% %% Load final data & fill with a complete trl_struct matrix.
% 
% % load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean_repairedchan/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal2'],'dataEEG_cleanedFinal2');
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean_repairedchan_keepchan/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal2'],'dataEEG_cleanedFinal2');
% 
% %make sure only EEG channels make it to the last dataset:
% cfg  = [];
% cfg.channel = {'all', '-Mast_L', '-Mast_R', '-HEOG_L','-HEOG_R', '-VEOG_Inf','-VEOG_Sup',  '-VEOG_Sp', '-ECG_Inf','-ECG_Sup'}; 
% dataEEG_cleanedFinal2 = ft_preprocessing(cfg, dataEEG_cleanedFinal2);
% 
% final_trialnums = dataEEG_cleanedFinal2.trialinfo(:,8);
% 
% final_trialrows = ismember(trl_struct(:,11), final_trialnums); 
% 
% final_trlstruct = trl_struct(final_trialrows,:);
% 
% %delete first 3 columns
% final_trlstruct2 = final_trlstruct; %first save original
% final_trlstruct(:,1:3) = []; %delete first 3 columns
% 
% dataEEG_cleanedFinal2.trialinfo = final_trlstruct;
% % dataEEG_cleanedFinal2.trialinfo2 = final_trlstruct2;
% 
% save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean_repairedchan2_keepchan/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal2'],'dataEEG_cleanedFinal2', '-v7.3');
% 
% end
% %% Test again by plotting ERPs
% 
% % temporalInfo_Artifacts.manual_postICA = dataEEG_cleanedFinal.cfg.previous.artfctdef;
% % temporalInfo_Artifacts.manual_postICA.channel = badchan;
% % save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/','sub',num2str(subjNumStr),'_SummaryRejectedArtifacts'],'temporalInfo_Artifacts');
% 
% %% Test plot 1: ERP
% clear;
% subjNumStr = '01';
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean_repairedchan2_keepchan/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal2'],'dataEEG_cleanedFinal2');
% % load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/trl_struct_files/sub',num2str(subjNumStr),'_trl_struct'],'trl_struct');
% 
% % load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/SemiautoClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedSemiauto'],'dataAll_EEGcleanedSemiauto');
% 
% datatemp = dataEEG_cleanedFinal2;
% 
% trl_struct(:,1:3) = [];
% final_trlstruct = trl_struct;
% % goodtrials = datatemp.cfg.trials;
% % datatemp = dataAll_EEGcleanedSemiauto;
% cfg = [];
% cfg.trl = final_trlstruct;
% % cfg.channel = {'Oz'};
% % cfg.channel = {'all', '-HEOG_L','-HEOG_R', '-VEOG_Inf','-VEOG_Sp','-VEOG_Sup', '-ECG_Inf','-ECG_Sup'}; %select only EEG channels
% cfg.channel = {'O1','Oz', 'Iz','O2','POO1','POO2','POz'};
% cfg.demean = 'yes';
% cfg.baselinewindow = [-0.5 0]; %min duration of prestim fixation point relativ to target image
% data_temp2 = ft_preprocessing(cfg, datatemp);
% 
% cfg.latency = [-0.5000 0.5000];
% dataEEG_cleanedFinal2 = ft_selectdata(cfg, data_temp2);
% 
% dataEEG_cleanedFinal2.trialinfo = final_trlstruct;
% 
% Normaltrials = find(dataEEG_cleanedFinal2.trialinfo(:,9) == 1);
% Maskedtrials = find(dataEEG_cleanedFinal2.trialinfo(:,9) == 2);
% LSFtrials = find(dataEEG_cleanedFinal2.trialinfo(:,9) == 3);
% 
% cfg = [];
% cfg.latency = [-.5 0.5];
% cfg.vartrllength = 2;
% cfg.trials = Normaltrials;
% test_normal = ft_timelockanalysis(cfg,dataEEG_cleanedFinal2);
% cfg.trials = Maskedtrials;
% test_masked = ft_timelockanalysis(cfg,dataEEG_cleanedFinal2);
% cfg.trials = LSFtrials;
% test_LSF = ft_timelockanalysis(cfg,dataEEG_cleanedFinal2);
% 
% figure;
% cfg = [];
% ft_singleplotER(cfg,test_normal,test_masked,test_LSF);
%         legend('Normal','Masked','LSF','Location','northwest')
% %         
% % %% Test plot 2: Time-frequency plot
% % 
% % cfg = [];
% % cfg.channel = {'O1','Oz', 'Iz','O2','POO1','POO2','POz'};
% % cfg.method = 'mtmfft';
% % cfg.output = 'pow';
% % cfg.foilim = [1 40];
% % cfg.tapsmofrq = 1;
% % cfg.taper = 'hanning';
% % cfg.vartrllength = 1;
% % cfg.trials = Normaltrials;
% % test_normal = ft_freqanalysis(cfg,dataEEG_cleanedFinal2);
% % cfg.trials = Maskedtrials;
% % test_masked = ft_freqanalysis(cfg,dataEEG_cleanedFinal2);
% % cfg.trials = LSFtrials;
% % test_LSF = ft_freqanalysis(cfg,dataEEG_cleanedFinal2);
% % 
% % figure;     
% % cfg = [];
% % ft_singleplotER(cfg,test_normal,test_masked,test_LSF);
% %     legend('Normal','Masked','LSF','Location','northeast')
