function results = CreateICA_Figures(dataEEG_ICA , numComps, subjNum, samplingFrequency)

cd '/isilon/LFMI/VMdrive/Ayaka/EEG/Preprocessing/';

icaComponentStruct = dataEEG_ICA;

results = C2F_createICAwindow(numComps, subjNum, samplingFrequency, icaComponentStruct);

path = pwd;

% folderName = ['ICAcomponents_',num2str(subjNum)] ;   % new folder name 
% myfolder=folderName;
% mkdir([path,filesep,myfolder]) ;
% path  = [path,filesep,myfolder] ;
% for c = 1:numComps
%     figure(c);
%     temp=[path,filesep,'Comp',num2str(c),'.png'];
%     saveas(gca,temp);
% end 

end
    