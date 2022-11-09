%% Preprocessing code [Load edf -> event import -> ASR -> epoch extraction -> baseline removal -> zscore ]
% 20-11-2017
% Madhavan
% first part automation: akshara
% modified: akshara 20 oct 2018
%% Adding eeglab to path
addpath('/home/aksharas/shared/softwares/eeglab14_1_1b');
eeglab;

%% User Inputs
phaseId = 2;
dataPath = sprintf('/home/data/EEG_Phase%d',phaseId);
% dataPath = sprintf('/home/data/EEG_Phase%d/pilotExpts',phaseId); %for pilot expts
savePath = sprintf('/home/aksharas/shared/dataResults/P%d/extracted_data',phaseId);

for subId = [1]
%% Load edf file 
fprintf('*************Loading data of Subject:%d************* \n',subId);
subName = sprintf('P%d_S%d',phaseId,subId);
%edfName = sprintf('%s/%s/filtered/%s.edf',dataPath,subName,subName);
% ^ use this if u have bandpass filtered data (band: 0.1- 70Hz; 50Hz notch)
% else do that later
edfName = sprintf('%s/%s/%s.edf',dataPath,subName,subName);
EEG = pop_fileio(edfName);
EEG.setname=subName; 
EEG = eeg_checkset( EEG );

%% event import
disp('*************Importing Events*************');
evtName = sprintf('%s/%s/%s.evt',dataPath,subName,subName);
EEG = pop_importevent( EEG, 'event',evtName,'fields',{'number' 'type' 'latency' 'urevent' 'duration'},'skipline',1,'timeunit',1);
% BE CAREFUL ABOUT THE NUMBER OF SKIPLINES (HEADER LINES TO SKIP)

eeglab draw % to update GUI to the current variable.
%% Load Channel locations
    disp('*************LOADING CHANNEL MAP*************');
    EEG=pop_chanedit(EEG, 'lookup','/home/aksharas/shared/softwares/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG,'nochannel',{'A1' 'A2' 'SECINDEX' 'SAMPINDEX'}); % removes last 4 channels
%     EEG = pop_select( EEG,'nochannel',{'A1' 'A2'});

eeglab draw % to update GUI to the current variable.
%% Downsample the data to 256 Hz 
% not doing for Phase 1 as we have not followed it in first analysis 
% Do it from phase 2 (Tools > Change Sampling rate)

%% filtering
% 1. band pass filter 0.1-70Hz 
% 1a. High pass filter with lower edge: 0.1Hz 
EEG = pop_eegfiltnew(EEG,  1);
% 1b. Low pass filter with higher edge: 70Hz
EEG = pop_eegfiltnew(EEG, [], 50);
%EEG = pop_eegfiltnew(EEG, [], 70, 194, 0, [], 1);
%pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,
%                                     revfilt, usefft, plotfreqz, minphase, 
%                                     usefftfilt);
% 2. notch filter at 50Hz 
% using CleanLine plugin doesnot seem to work always.Parameters:step size==window size(def. 4s) 
EEG = pop_eegfiltnew(EEG, 49.5, 50.5, [], 1);
%EEG = pop_eegfiltnew(EEG, 49.5, 50.5, 6760, 1, [], 1);

%% Remove bad channels (P1: used prep; makato (eeglab) suggests clean_rawdata-> P2)
% [nCh,nSamples] = size(EEG.data);
% load('/home/aksharas/shared/dataResults/P1R1/extracted_data/bad_channels.mat');
% badCh = bad_channels{subId};
% 
% for i=1:length(bad_channels{1,subId})
%     x=bad_channels{1,subId};
%     z=char(x(1,i));
%     EEG = pop_select( EEG,'nochannel',{z}); 
% end

% Remove breaks

% %% ASR
EEG = clean_rawdata(EEG, -1, [-1],-1, -1, 12, 'off');
%^ without channel removal
% EEG = clean_rawdata(EEG, 'off', [0.25 0.75], 'off', 'off', 5, 0.5);
% EEG = eeg_checkset( EEG );

%% Epoch extraction
saveFold =  sprintf('%s/%s',savePath,subName);
if (~exist(saveFold,'dir'))
    mkdir(saveFold);
end

% Extract rest state data
restEEG = pop_epoch( EEG, {'rest'}, [0 1.5], 'newname', sprintf('%s_rest',subName), 'epochinfo', 'yes');
saveName = sprintf('%s_rest.set',subName);
restEEG = pop_saveset( restEEG, 'filename',saveName,'filepath',saveFold);

% Extract prestim period data
baselineEEG = pop_epoch( EEG, {'prestim'}, [0 0.4], 'newname', sprintf('%s_baseline',subName), 'epochinfo', 'yes');
saveName = sprintf('%s_baseline.set',subName);
baselineEEG = pop_saveset( baselineEEG, 'filename',saveName,'filepath',saveFold);

% Extract listening state data
listeningEEG = pop_epoch( EEG, {'listening'}, [0 1.5], 'newname', sprintf('%s_listening',subName), 'epochinfo', 'yes');
saveName = sprintf('%s_listening.set',subName);
listeningEEG = pop_saveset( listeningEEG, 'filename',saveName,'filepath',saveFold);

% Extract speaking state data
speakingEEG = pop_epoch( EEG, {'speaking'}, [0 1.5], 'newname', sprintf('%s_speaking',subName), 'epochinfo', 'yes');
saveName = sprintf('%s_speaking.set',subName);
speakingEEG = pop_saveset( speakingEEG, 'filename',saveName,'filepath',saveFold);

%% Remove bad trials of listening state EEG
%{
[numChannels,numSamples,numEpochs] = size(listeningEEG.data);
load('/home/aksharas/shared/dataResults/extracted_data/bad_trials.mat');
% load(sprintf('%s/bad_trials.mat',savePath));  %numSub x 1 cell containing index of bad trials (for 560x1) %CHANGE loc
%[should remove before specific stimulus as index corresponds to 560x1 
disp('Removing bad trials......');
[~,~,numEpochs] = size(listeningEEG.data); 
badInd = ones(numEpochs,1);

badInd(bad_trials{subId}') = 0;
data = listeningEEG.data(:,:,logical(badInd));
[~,~,numEpochs] = size(data); %update size 

load('/home/aksharas/shared/dataResults/extracted_data/prompts_wordsPhrases.mat');
prompts = prompts(logical(badInd)); %remove bad trials from prompts too.
% save(sprintf('%s/prompts_wordsPhrases_badTrialsRemoved.mat',saveFold),'prompts');
%}

%% Baseline removal
%Compute baseline average
baseline_avg = (mean(baselineEEG.data,2));

load(sprintf('%s/prompts.mat',savePath));

%baseline_avg_wp = baseline_avg(:,:,logical(badInd));

%remove baseline average from each state data
[numChannels,numSamples,numEpochs] = size(listeningEEG.data);
 data = listeningEEG.data;
 
listenData = data-baseline_avg;
% imagineData = imaginingEEG.data-baseline_avg_wp; % will do later after
% finding bad_trials
% speakData = speakingEEG.data-baseline_avg_wp;
% 
% %remove baseline average from each state data of sentences
% sent_listenData = sent_listeningEEG.data-baseline_avg_s;
% sent_imagineData = sent_imaginingEEG.data-baseline_avg_s;
% sent_speakData = sent_speakingEEG.data-baseline_avg_s;

%% Z-score computation
    listenData = reshape(listenData,[numChannels, numSamples*numEpochs]);
    listenData = zscore(listenData);
    listenData = reshape(listenData,[numChannels, numSamples,numEpochs]);
%load('/home/aksharas/shared/dataResults/P1R1/extracted_data/bad_channels_number.mat');
%% Introduce NaNs back into channels 
    bad_ch=bad_channels_number{1,subId};
for i=1: length(bad_channels_number{1,subId})
    row_no=bad_ch(1,i); %%where wants to insert
    listenData(1:row_no-1,:,:) = listenData(1:row_no-1,:,:);
    tp =listenData(row_no:end,:,:);
    listenData(row_no,:,:)=NaN;
    if (row_no~=30)
    listenData(row_no+1:end+1,:,:) =tp;
    end
end
%% Introduce NaNs back into trials
  bad_tr=bad_trials{1,subId};
for i=1: length(bad_trials{1,subId})  
    row_no=bad_tr(1,i); %%where wants to insert
    listenData(:,:,1:row_no-1) = listenData(:,:,1:row_no-1);
    tp =listenData(:,:,row_no:end);
    listenData(:,:,row_no)=NaN;
    if (row_no~=length(prompts)) %
    listenData(:,:,row_no+1:end+1) =tp;
    end
end

%% 
save(sprintf('%s/%s_listenZscore.mat',saveFold,subName),'listenData');

end