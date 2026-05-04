function plotElectrodePowerSpectrums(data,labels,channelIndices,fs,subject)
%% Parameters
% data= EEGdata.trial{1}; 
% labels= EEGdata.label;
% channelIndices= 1:nchans; %, powerSpectrumDirectory, EEGdata.fsample,'_RAW'); % could add topoplot to figures here
% fs=EEGdata.fsample;

nFFT = 2^12;
numChans = size(data,1);
noFigs = ceil(numChans/9);
psds = struct;

[p,f] = pwelch(data(1,:),hamming(2^12),[],2^12,fs);
ps = zeros(numChans,length(p));


for i_chan = 1:numChans
    ps(i_chan,:) =  pwelch(data(i_chan,:),hamming(nFFT),[],nFFT,fs);
end

averagePowerSpectrum = squeeze(mean(log10(ps(find(ps(:,1)~=0),:))));

% for fg = 1:noFigs
%     figure
%     k = 1;
%     for i_chan = [1 + (fg-1)*9:min(numChans,fg*9)]
%         h(k) = subplot(3,3,k);
%         plot(f,squeeze(log10(ps(i_chan,:))));
%         hold on
%         plot(f,averagePowerSpectrum,'linewidth',2,'color','k');
%         k= k+1;
%         title(['Chan ' num2str(channelIndices(i_chan)) ' ' labels{i_chan}]);
%         set(gca,'xlim',[0 120]);
%     end
%     %set(gcf,'position',[ 302          85        1664         981]);
%     %print('-dpng',[outputFigureDirectory 'PowerSpectrum' num2str(fg) suffix '.png']);
%     %close all
% end

for fg = 1:noFigs
    figure; %('visible', 'off');
    k = 1;
    for i_chan = [1 + (fg-1)*9:min(numChans,fg*9)]
        h(k) = subplot(3,3,k);
        plot(f,squeeze(log10(ps(i_chan,:))));
        hold on
        plot(f,averagePowerSpectrum,'linewidth',2,'color','k');
        k= k+1;
        title(['Chan ' num2str(channelIndices(i_chan)) ' ' labels{i_chan}]);
        set(gca,'xlim',[0 50]);
    end
%     set(gcf,'position',[ 302          85        1664         981]);
%     print('-dpng',[outputFigureDirectory 'LowFreqPowerSpectrum' num2str(fg) suffix '.png']);
%     close all
end
path = pwd ;
folderName = ['ChannelPwrSpct_',num2str(subject)] ;   % new folder name 
myfolder=folderName;
folder = mkdir([path,filesep,myfolder]) ;
path  = [path,filesep,myfolder] ;
for c = 1:noFigs
    figure(c);
    plot(rand(1,10));
    temp=[path,filesep,'fig',num2str(c),'.png'];
    saveas(gca,temp);
end
