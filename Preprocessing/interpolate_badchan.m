%% Interpolate bad channels.
%(Repair channels)
%preparing layout

%% This version of the "inerpolate bad channels" tries to delete as few channels as possible, unless it was very obviously dead during the entire recording.
% Theoretically, this should improve MVPA.

clear; clc

paths_C2F_EEG = C2F_EEG_SetPaths;
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
vars = vars_C2F_EEG;
addpath('/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing');

for subjNum = [39] %[1 2 3 4 7 8 10 11 13 14 15 17 18 20 22 24 25 26 27 28]
subjNumStr = SubLabel_Add0BeforeSingleDigit(subjNum); % add a 0 before single digit subjNums

if subjNum == 18 || subjNum == 27
    origdata = load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FilteredOnly/','sub',num2str(subjNumStr),'a_dataAll_EEGfiltered'],'data_p3');
else
    origdata = load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FilteredOnly/','sub',num2str(subjNumStr),'_dataAll_EEGfiltered'],'data_p3');
end

finaldata = load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/ICAclean/','sub',num2str(subjNumStr),'_dataEEG_cleanedICA'],'dataEEG_cleanedICA');


if subjNum == 33 || subjNum == 38
    cfg.channel = {'all', '-Mast_L', '-Mast_R', '-HEOG','-VEOG', '-ECG'}; %only EEG channels
    cfg.reref = 'yes';
    cfg.refchannel = {'all'};
    cfg.refmethod = 'avg';
    finaldata.dataEEG_cleanedICA = ft_preprocessing(cfg,finaldata.dataEEG_cleanedICA);
end

origchan = origdata.data_p3.label;
finalchan = finaldata.dataEEG_cleanedICA.label;
bad_chan = setdiff(origchan, finalchan);

badchan_all = [];
for c = 1:length(bad_chan)
    temp = bad_chan{c};
    if contains(temp, 'Mast_L') || contains(temp, 'Mast_R') || contains(temp, 'HEOG_L') || contains(temp, 'HEOG_R') || contains(temp, 'VEOG_Inf') || contains(temp, 'VEOG_Sp')
        badchan_all = badchan_all;
    else
        badchan_all{c} = temp;
    end
end

if ~isempty(badchan_all)
    badchan_all = badchan_all(~cellfun('isempty',badchan_all));
end
 
load 'acticap128.mat'   
elec = struct;
for i = 3:length(acticap128.chanlocs)% excluding the ground and reference (1+2)
elec.pnt(i-2,1) = acticap128.chanlocs(i).X;
elec.pnt(i-2,2) = acticap128.chanlocs(i).Y;
elec.pnt(i-2,3) = acticap128.chanlocs(i).Z;
elec.label{i-2} = acticap128.chanlocs(i).labels;
end

cfg = [];
cfg.elec = elec;
cfg.rotate = 90; 
cfg.layout     =ft_prepare_layout(cfg);
cfg.template      = 'easycap128.mat';
cfg.method = 'triangulation';
cfg.neighbours = ft_prepare_neighbours(cfg);   
neighbours= cfg.neighbours;

cfg.badchannel = badchan_all; % {''};
% cfg.badchannel = {'Fp2', 'F3', 'O1', 'AF3', 'FPz', 'FC6', 'Fz', 'FFC4h'}; 
% cfg.badchannel = {'FC5', 'CP2', 'TP10', 'FC4', 'TP8', 'FCC4h', 'FTT10h', 'CPP3h', 'POO9h', 'FFT9h', 'FFT10h'};
% cfg.missingchannel = sO9','Ol1h', 'Ol2h', 'O10', 'P9', 'P10'};% these channels were missing for all participants
cfg.method = 'average';
cfg.neighbours = neighbours;
cfg.trials = 'all';
%cfg.lambda = regularisation parameter (default = 1e-5, only for methods "spline" and "slap")
%cfg.order = 4; % (default, only for methods "spline" and "slap")
dataEEG_cleanedFinal2 = ft_channelrepair(cfg, finaldata.dataEEG_cleanedICA);

save(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/FinalClean_repairedchan_keepchan/','sub',num2str(subjNumStr),'_dataEEG_cleanedFinal2'],'dataEEG_cleanedFinal2', '-v7.3');

end