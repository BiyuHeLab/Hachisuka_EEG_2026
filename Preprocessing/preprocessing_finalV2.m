%% Preprocess EEG data
% This script includes all steps of preprocessing for EEG data
    %1. load in subject and task wise EEG/EOG/ECG-data
    %2. Determine trials (or blocks) from the continious EEG recording by means of trigger file
    %3. Apply basic preprocessing and re-referencing to EEG data
    %4. Apply basic preprocessing and re-referencing to EOG & ECG data
    %5. Append EEG and EOG/ECG data to one file
    %6. Semi-automatically check for jump- and muscle artifacts
    %7. Manual artifact rejection 1 - determine broken/noisy channels and remaining jumps
    %8. ICA & MI (mutual information) approach - compute an ICA on the EEG
        %data and then check mutual information with separate EOG/ECG channels
        %in order to remove right components, then remove said components
        %and backproject to channel level
    %9. Manual artifact rejection 2 - check data quality after ICA/MI and
        %remove any remaining noisy trials/channels
    %10. Repair any removed channels by extrapolating from adjecent
        %channels
        
%NOTES:
    %8 Missing electrodes (used as EOG/ECG): 
    %PPO2h, PPO6h, PPO10h, P10, I1, Ol1h, Ol2h, I2, PPO2h, PPO6h, PPO10h,P10
    %Previously listed but incorrect: O9,Ol1h, Ol2h, O10, P9, P10
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

paths_C2F_EEG = C2F_EEG_SetPaths;
vars = C2F_EEG_SetVars(paths_C2F_EEG);
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
cd /isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Scripts
initDir = pwd;

ft_default.trackusage = 'no';

% subjNum = vars.validSubjs; %Specify current subject
subjNum = 37;
subjNumStr = SubLabel_Add0BeforeSingleDigit(subjNum); % add a 0 before single digit subjNums

%1. Load in behavioral data
% load([vars.C2F_behaviorDir, 'C2F_behavData.mat']);
rawBehavData = C2F_EEG_createBehavData(subjNum);
% Load in EEG data
taskTypeNum = 2; %because only 2(Experiment) has EEG data
taskType = vars.Task.TaskTypes{taskTypeNum};

%Specify EEG data path
subjDir = ['sub' subjNumStr '/Experiment_data/' taskType  '/'];
% rawDataFileName = ['C2F_EEG_00' subjNumStr '_2.vhdr'] ; % ['sub' subjNumStr 'pos_allBlocks.vhdr'];
rawDataFileName = ['C2F_EEG_00' subjNumStr '.vhdr'] ; % ['sub' subjNumStr 'pos_allBlocks.vhdr'];
preProcessingSettingsName = 'default_EEG';
eegDataDir = ['/isilon/LFMI/VMdrive/Ayaka/EEG/fromThomas/Data/RawData/', subjDir,'EEG_rawData/'];
cd(eegDataDir)
RawEEGDataFileName = [eegDataDir rawDataFileName];    
EEGTriggerFileName = [rawDataFileName(1:end-4) 'vmrk'];

% Preprocessing begins here:
cfg = [];
cfg.dataset = RawEEGDataFileName;
cfg.continuous  = 'yes'; %Tell FT that the data is one continuous recording
% cfg.trialfun = 'C2F_trialfun_defineblocks';
cfg.trialfun = 'C2F_trialfun_definetrials';
cfg = ft_definetrial(cfg);
cfg.hdr = ft_read_header(cfg.dataset);
%%
%2.1 Read out trial-specific (i.e., condition) and behavioral (i.e., response) data and append it to the trialinfo subfield,
%determine fit by checking/comparing the start times of the respective trial in the EEG and behavioral file
rawBehavData2(subjNum) = rawBehavData;
rawBehavData = rawBehavData2;
for block = 1:length(rawBehavData(subjNum).Experiment.blockData)
    for trial_perblock = 1:length(rawBehavData(subjNum).Experiment.blockData(block).trialData.trialnumber)
        
        trial_starttime = rawBehavData(subjNum).Experiment.blockData(block).triggerData.StartTime_block...
                                   +rawBehavData(subjNum).Experiment.blockData(block).triggerData.StartTimeAligned2BlockStart_startphase(trial_perblock);
        trial = find(cfg.trl(:,1) == trial_starttime);        
                   
                %Condition               
                if strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.condition(trial_perblock),'Normal')                      
                cfg.trl(trial,12) = 1;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.condition(trial_perblock),'Masked')   
                cfg.trl(trial,12) = 2;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.condition(trial_perblock),'LSF') 
                cfg.trl(trial,12) = 3;
                end
                
                %Superordinate
                if strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.domain(trial_perblock),'Animate')                      
                cfg.trl(trial,13) = 1;
                elseif strcmp(rawBehavData(subjNum).Experiment.blockData(block).trialData.domain(trial_perblock),'Inanimate')   
                cfg.trl(trial,13) = 2;
                end 
                
                %Basic                
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

trl_struct = cfg.trl; %store trialstructure in separate variable for later.
save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/trl_struct_files/sub',num2str(subjNumStr),'_trl_struct'],'trl_struct', '-v7.3');

% cfg.continuous = 'no'; %to view raw data by blocks.
% ft_databrowser(cfg);
%% 2. Preprocessing
% % This step will read, filter and cut the raw data
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/trl_struct_files/sub',num2str(subjNumStr),'_trl_struct'],'trl_struct');
% Because we are rereferencing to mastoids, first check their quality.
% If mastoid recording problematic, then re-referencing with them doesn't
% make sense!

block_col = cfg.trl(:,8);
% for block = (1:length(rawBehavData(subjNum).Experiment.blockData)) 
for block = 1:block_col(end)
        cfg2 = [];
        cfg2.dataset = RawEEGDataFileName;
        cfg2.continuous  = 'yes'; %Tell FT that the data is one continuous recording
        cfg2.trl = cfg.trl([cfg.trl(:,8) == block],:);
        cfg2.channel = {'Mast_L', 'Mast_R'};    
        test = ft_preprocessing(cfg2);
        cfg2 = [];
        cfg2.method = 'channel';
        cfg2.keeptrial = 'yes';
        cfg2.keepchannel = 'yes'; %Keep noisy channels. wel'll delete them later
        ft_rejectvisual(cfg2,test);  
end
clear test   
                
% 2.1. Read in the raw data and re-reference EOG & EEG:

cfg.channel = {'all'}; %EOG & EEG 
cfg.reref = 'yes';
if subjNum == 33
    cfg.reref = 'no';
    data_p = ft_preprocessing(cfg);
else
    if subjNum == 17 || subjNum == 36 
        cfg.refchannel = {'Mast_R'};
    else
        cfg.refchannel = {'Mast_L', 'Mast_R'};
    end
    cfg.refmethod = 'avg';
    data_p = ft_preprocessing(cfg);
end

% store EOG channels to append later:
cfg = [];
cfg.continuous  = 'yes';
cfg.trl = trl_struct; %Use previously specified trial structure
cfg.channel = {'HEOG_L','HEOG_R', 'VEOG_Inf','VEOG_Sup',  'VEOG_Sp'};
cfg.continuous  = 'yes'; %preprocess the entire recording.
cfg.demean = 'yes'; % Remove mean across time to get rid of offset.
cfg.detrend = 'yes'; %Fit a line and subtracting the fit.
cfg.hpfilter = 'yes'; % Enable high-pass filter to remove slow drifts
cfg.lpfilter = 'yes';  % Enable low-pass filter
cfg.hpfiltord = 3;
cfg.lpfiltord = 3;
cfg.hpfilttype = 'but';
cfg.lpfilttype = 'but';
cfg.hpfreq = 0.05; % set up the frequency for high-pass filter
cfg.lpfreq = 150; % set up the frequency for low-pass filter
cfg.bsfilter = 'yes';    
cfg.bsfilttype = 'but'; %Bandstop filter for line noise
cfg.bsfreq = [59 61 119 121];
dataEOG = ft_preprocessing(cfg, data_p);
 
% cfg =[];
% cfg.method = 'trial';
% data_EOG = ft_rejectvisual(cfg, dataEOG);

% cfg = [];
% cfg.continuous = 'no'; %to view by blocks.
% ft_databrowser(cfg, data_p);

cfg = [];
cfg.channel = {'all', '-HEOG_L','-HEOG_R', '-VEOG_Inf','-VEOG_Sup',  '-VEOG_Sp', '-ECG_Inf','-ECG_Sup', '-Mast_L', '-Mast_R'};  %EEG only
cfg.continuous  = 'yes'; %preprocess the entire recording.
cfg.demean = 'yes'; % Remove mean across time to get rid of offset.
cfg.detrend = 'yes'; %Fit a line and subtracting the fit.
cfg.hpfilter = 'yes'; % Enable high-pass filter to remove slow drifts
cfg.lpfilter = 'yes';  % Enable low-pass filter
cfg.hpfiltord = 3;
cfg.lpfiltord = 3;
cfg.hpfilttype = 'but';
cfg.lpfilttype = 'but';
cfg.hpfreq = 0.05; % set up the frequency for high-pass filter
cfg.lpfreq = 150; % set up the frequency for low-pass filter
cfg.bsfilter = 'yes';    
cfg.bsfilttype = 'but'; %Bandstop filter for line noise
cfg.bsfreq = [59 61 119 121];
data_p2 = ft_preprocessing(cfg, data_p);

% cfg = [];
% cfg.continuous = 'no'; %to view by blocks.
% cfg.method = 'trial';
% ft_rejectvisual(cfg, data_p2);

%% Add trial_info header
data_p3 = ft_appenddata([], data_p2, dataEOG);
data_p3.trialinfo_header =  {'responseslide_superordinate_offsettotrialstart', 'responseslide_basic_offsettotrialstart'...
                                    ,'responseslide_PAS_offsettotrialstart', 'trial_eyetrackingsample'...
                                    ,'block_number', 'block_starttime', 'trial_numberinblock','trial_totalnumber'...
                                    ,'condition (1 = Normal, 2 = Masked, 3 = LSF)','domain (1 = Animate, 2 = Inanimate)'...
                                    ,'class (1 = Dog, 2 = Cat, 3 = Car, 4 = Truck'...
                                    ,'targetimage (1:3 = Dog1:3, 4:6 = Cat1:3, 7:9 = Car1:3, 10:12 = Truck1:3)','contrast_level','image_location (1 = Top, 2 = Bottom, 3 = Left, 4 = Right)'...
                                    ,'responsecorrectness_superordinate (0 = false, 1 = correct)','RT_Superordinate','resposnecorrectness_basic (0 = false, 1 = correct)','RT_basic'...
                                    ,'response_PAS (0:3)','RT_PAS'...
                                    ,'duration_startphase','duration_prestimjitter','duration_targetimage','duration_delaymask','duration_poststimjitter'};       

% save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FilteredOnly/','sub',num2str(subjNumStr),'_dataAll_appended'],'data_p3', '-v7.3');
save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FilteredOnly/','sub',num2str(subjNumStr),'_dataAll_EEGfiltered'],'data_p3', '-v7.3');

%% 6. Preprocessing Step 1 - remove jump and muscle artifacts (important to do before ICA)
%Aim; Remove big/prominent jump and muscle artifacts - be liberal, since we will do the thorough rejection afterwards
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FilteredOnly/','sub',num2str(subjNumStr),'_dataAll_EEGfiltered'],'data_p3');

%If it's a particularly bad subject, manually reject channels here:
% cfg = []; % empty the cfg-structure
% cfg.method = 'trial';
% cfg.channel =  {'all', '-Mast_L', '-Mast_R'};
% data_p3 = ft_rejectvisual(cfg,data_p3); % we call the databrowser with the cfg-settings and the dat-structure defined above

%6.1Jumps
cfg = [];
cfg.trl = trl_struct; %separate into trials, because we want to remove trials when there are huge artifacts.
cfg.continuous = 'no'; %because we have trial-wise data now

cfg.artfctdef.jump.channel = {'all', '-Mast_L', '-Mast_R', '-HEOG_L','-HEOG_R', '-VEOG_Inf','-VEOG_Sup',  '-VEOG_Sp', '-ECG_Inf','-ECG_Sup'}; %only EEG channels
cfg.artfctdef.jump.cutoff = 150;
cfg.artfctdef.jump.trlpadding = 0;
cfg.artfctdef.jump.artpadding = 0;
cfg.artfctdef.jump.fltpadding = 0;
cfg.artfctdef.jump.cumulative = 'yes';
cfg.artfctdef.jump.medianfilter = 'yes';
cfg.artfctdef.jump.medianfilterord = 9;
cfg.artfctdef.jump.absdiff = 'yes';
cfg.artfctdef.jump.interactive = 'yes';

[cfg, artifcat_jump] = ft_artifact_jump(cfg, data_p3);

%6.2Muscle
cfg.artfctdef.muscle.channel = {'all', '-Mast_L', '-Mast_R', '-HEOG_L','-HEOG_R', '-VEOG_Inf','-VEOG_Sup',  '-VEOG_Sp', '-ECG_Inf','-ECG_Sup'};  %only EEG channels
cfg.artfctdef.muscle.cutoff = 80;
cfg.artfctdef.muscle.trlpadding = 0;
cfg.artfctdef.muscle.artpadding = 0;
cfg.artfctdef.muscle.fltpadding = 0.1;
cfg.artfctdef.muscle.bpfilter = 'yes';
cfg.artfctdef.muscle.bpfreq = [110 140];
cfg.artfctdef.muscle.bpfiltord = 9;
cfg.artfctdef.muscle.bpfilttype = 'but';
cfg.artfctdef.muscle.hilbert = 'yes';
cfg.artfctdef.muscle.boxcar = 0.2;
cfg.artfctdef.muscle.interactive = 'yes';

[cfg, artifcat_muscle] = ft_artifact_muscle(cfg, data_p3);

%6.3 Remove semi-automatically determined artifacts
% cfg.artfctdef.feedback = 'yes';
cfg.artfctdef.reject = 'complete';
% cfg.artfctdef.minaccepttim = 2.2; %Minimum length in s of remaining trial (2.2s without response)
dataAll_EEGcleanedSemiauto = ft_rejectartifact(cfg,data_p3);

save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/SemiautoClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedSemiauto'],'dataAll_EEGcleanedSemiauto', '-v7.3');

%% Check data again & remove still-broken channels after filtering & artifact removal:
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/SemiautoClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedSemiauto'],'dataAll_EEGcleanedSemiauto');

%Note: Prior to ICA delete those bad channels that are dead/flat or corrupted, 
%but NOT those that have those artifacts that I want to clean with ICA
%(i.e. cardiac, EOG).

cfg = []; % empty the cfg-structure
cfg.method = 'summary';
% cfg.channel = {'all', '-Mast_L','-Mast_R','-ECG'}; % EEG & EOG channels only, remove mastoid channels.
cfg.channel =  {'all', '-Mast_L', '-Mast_R'};
% cfg.channel = {'all', '-Mast_L', '-Mast_R', '-HEOG','-VEOG', '-ECG'}; %only EEG channels
dataAll_EEGcleanedManually = ft_rejectvisual(cfg,dataAll_EEGcleanedSemiauto); % we call the databrowser with the cfg-settings and the dat-structure defined above
bad_chan1 = setdiff(dataAll_EEGcleanedSemiauto.label, dataAll_EEGcleanedManually.label);

save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/ManualClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedManually'],'dataAll_EEGcleanedManually', '-v7.3');

%save badchan_part1 and bad_chanpart2 for post-ICA faulty electrode
%interpolation:
% badchans = bad_chan2;
badchans = unique([bad_chan1]);
save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/ManualClean/badchans/','sub',num2str(subjNumStr),'_badchans'],'badchans', '-v7.3');

disp(badchans);

% optional: use "var" to exclude excessively bad trials:
cfg = []; % empty the cfg-structure
cfg.method = 'trial';
% cfg.channel = {'all', '-Mast_L','-Mast_R','-ECG'}; % EEG & EOG channels only, remove mastoid channels.
cfg.channel =  {'all', '-Mast_L', '-Mast_R'};
dataAll_EEGcleanedManually = ft_rejectvisual(cfg,dataAll_EEGcleanedManually); % we call the databrowser with the cfg-settings and the dat-structure defined above

save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/ManualClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedManually'],'dataAll_EEGcleanedManually', '-v7.3');

bad_chan2 = setdiff(dataAll_EEGcleanedSemiauto.label, dataAll_EEGcleanedManually.label);
badchans = unique([bad_chan2]);
save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/ManualClean/badchans/','sub',num2str(subjNumStr),'_badchans'],'badchans', '-v7.3');

disp(badchans);

% IF SUBJ 33, rereference to average:
% cfg.channel = {'all'}; %EOG & EEG 
% cfg.reref = 'yes';
% cfg.refchannel = {'all'};
% cfg.refmethod = 'avg';
% dataAll_EEGcleanedManually = ft_preprocessing(cfg,dataAll_EEGcleanedManually);

%% Aggressively filter before ICA:
% load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/ManualClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedManually'],'dataAll_EEGcleanedManually');

cfg = [];
cfg.channel = {'all', '-3ast_L', '-Mast_R'};  %EEG & EOG only
cfg.continuous  = 'yes'; %preprocess the entire recording.
cfg.hpfilter = 'yes'; % Enable high-pass filter to remove slow drifts
cfg.hpfiltord = 3;
cfg.hpfilttype = 'but';
cfg.hpfreq = 1;
data_processedforICA = ft_preprocessing(cfg, dataAll_EEGcleanedManually);

save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'ICA/','sub',num2str(subjNumStr),'_data_processedforICA'],'data_processedforICA', '-v7.3');

% We can now stitch the segmented data back together in a continuous
% representation -- no need to do this, ICA on epoched data is preferred
% (will concatenate anyway), as long as there are enough datapoints.

% cfg = [];
% cfg.continuous = 'yes';
% dataAll_EEGcleanedManually_cont = ft_redefinetrial(cfg, dataAll_EEGcleanedManually);
% save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/ManualClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedManually_cont'],'dataAll_EEGcleanedManually_cont', '-v7.3');

