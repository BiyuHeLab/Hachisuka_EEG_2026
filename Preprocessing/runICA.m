%Test run for Big Purple:
%Add paths
clear;
clc;

cd '/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing'

paths_C2F_EEG = C2F_EEG_SetPaths;
vars = C2F_EEG_SetVars(paths_C2F_EEG);
vars_C2F_EEG = C2F_EEG_SetVars(paths_C2F_EEG);
initDir = pwd;

subjNum = 39;
subjNumStr = '39';

%load rawBehavData
rawBehavData = C2F_EEG_createBehavData(subjNum);

% subjNumStr = SubLabel_Add0BeforeSingleDigit(subjNum); % add a 0 before single digit subjNums

load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'Clean/ManualClean/','sub',num2str(subjNumStr),'_dataAll_EEGcleanedManually'],'dataAll_EEGcleanedManually');
load(['/isilon/LFMI/VMdrive/Ayaka/EEG/Data/EEG_preproc/', 'ICA/','sub',num2str(subjNumStr),'_data_processedforICA'],'data_processedforICA');

cd '/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing'

data_orig = dataAll_EEGcleanedManually;
data_aggressive = data_processedforICA;

% C2F_ICA_bigpurple(subjNum,subjNumStr, rawBehavData,dataAll_EEGcleanedManually);
C2F_ICA_bigpurple(subjNum, subjNumStr, rawBehavData, data_orig, data_aggressive)