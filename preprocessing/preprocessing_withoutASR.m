clc; clear all; close all;
%% User Inputs
phaseId = 2;
% parentPath = '/Volumes/AKSHARA'; %for pendrive
parentPath = '/home/user/Documents/projects'; %for labSys


dataPath = sprintf('%s/semantics/data/main',parentPath);
% dataPath = sprintf('/home/data/EEG_Phase%d/pilotExpts',phaseId); %for pilot expts
% savePath = sprintf('/Volumes/AKSHARA/semantics/dataResults/P%d/extracted_data',phaseId);
savePath = [parentPath '/semantics/preprocessing/temp_IS'];
mkdir(savePath);
chIndx = 1:56;
lf = 0.1; hf = 30; %lower and higher edges of passband
%% Adding eeglab to path
addpath('~/Documents/softwares/eeglab14_1_2b');
eeglab;

for subId = 1:20
%% load data
fprintf('*************Loading data of Subject:%d************* \n',subId);
subName = sprintf('P%d_S%d',phaseId,subId);
edfName = sprintf('%s/%s/%s.edf',dataPath,subName,subName);
EEG = pop_fileio(edfName);
EEG.setname=subName; 
EEG = eeg_checkset( EEG );

% Load Channel locations
    disp('*************LOADING CHANNEL MAP*************');
    EEG=pop_chanedit(EEG, 'lookup','~/Documents/softwares/eeglab14_1_2b/plugins/dipfit3.0/standard_BESA/standard-10-5-cap385.elp');
    EEG = eeg_checkset( EEG );
    EEG = pop_select( EEG,'nochannel',{'A1' 'A2', 'P2','Pz','PO4','PO8'}); % removes non-mapped channels


%% Step 1: Notch filter at 50Hz using narrow band-feedback structure
load('notchFilter/FBnotchCoeffs.mat');
datay = filter(num,den,EEG.data');
EEG.data = datay';

%% Step 2: FIR HPF at 0.1Hz
% Using IIR Butterworth filter (ERPLAB) of order=2 (mild)
EEG  = pop_basicfilter( EEG,  chIndx , 'Boundary', 'boundary', 'Cutoff',  lf, 'Design', 'butter', 'Filter', 'highpass', 'Order',  2, 'RemoveDC',...
 'on' );

%% Step 3: FIR LPF at 30Hz
% Using IIR Butterworth filter (ERPLAB) of order=4 (mild) OR 20Hz & order=2
EEG  = pop_basicfilter( EEG,  chIndx , 'Boundary', 'boundary', 'Cutoff',  hf, 'Design', 'butter', 'Filter', 'lowpass', 'Order',...
  4 ); 

filename = sprintf('P2_S%d_afterFilters.set',subId); 
EEG = pop_saveset( EEG, 'filename',filename,'filepath',savePath);
% figure; pop_spectopo(EEG, 1, [0      5520999.0234], 'EEG' , 'percent', 15, 'freq', [6 10 22], 'freqrange',[0 120],'electrodes','off');
% %% load pre-processed data 
% EEG = pop_loadset('filename',sprintf('P2_S%d_afterFilters.set',subId),'filepath',dataPath);

originalEEG = EEG; % to retreive org channel locs

    
%% Step 5: Eyeblink removal : ASR (clean_rawdata only for ASR)
 EEG = clean_rawdata(EEG, -1, [-1], 0.8, -1, 5, 0.5);  
filename = sprintf('P2_S%d_afterASR.set',subId); 
EEG = pop_saveset( EEG, 'filename',filename,'filepath',savePath);

% EEG = pop_loadset('filename',sprintf('P2_S%d_afterASR.set',subId),'filepath',dataPath);
if exist('EEG.etc.clean_channel_mask')
bad_channels_ASR = find(EEG.etc.clean_channel_mask==0); %to know the channels rejected by ASR
end     

%% Step 4: Remove Bad channels : Manually 
% 1. By time domain 2. Spectrum 
%>> Save this file to _afterASR.set

%% Step 6: Interpolate Removed Channels
EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
filename = sprintf('P2_S%d_afterInterp.set',subId); 
EEG = pop_saveset( EEG, 'filename',filename,'filepath',savePath);

%% Step 7: Re-referencing & Z-score computation(?) > will do later SAVE till step 6


%% Step 8: Bad trial removal : automatic and manual check


%% Epoching


%% Detrending ('ALL')

end


