%plots threshold vs sample for all the range time verification algorithms

close all;
clear all;

% Input variables 
Length = 500000; %size of the data (how many samples of noise)
PFA = 10^-2;    %Probability of false alarm
RefWindow = 32; %total window size (is divided in 2 for leading and lagging)
guardCells = 2; %total number of guard cells (is divided in 2 for leading and lagging)(must be even number greater than 0)
aGO = GOCA_CFAR_Alpha(PFA, RefWindow);
aSO = SOCA_CFAR_Alpha(PFA, RefWindow);

%Edgeclutter?
edgeClutter = 1;

%Noise Generator
samplesignal = zeros([1  Length]); %signal initialisation
signal = samplesignal;
t = 1:1:length(samplesignal); %time of simulation
samplesignal = normrnd(0,1,1,Length) + 1i*normrnd(0,1,1,Length);  %complex white noise generation

if edgeClutter == 1
    signal(1:200) = samplesignal(1:200);
    signal(200:end) = normrnd(2,1,1,Length-199) + 1i*normrnd(2,1,1,Length-199);
else
    signal = samplesignal;
end

DataAfterPowerLawDetector = abs(signal).^2; %realising signal power
if edgeClutter == 2
    DataAfterPowerLawDetector(1:200) = DataAfterPowerLawDetector(1:200);
    DataAfterPowerLawDetector(200:end) = DataAfterPowerLawDetector(200:end) +10;
end



%CA_CFAR
T_CA_CFAR  = zeros([1  length(signal)]);  %initialise an array for threshold values

for CUT = 1: length(signal)
    if CUT <= RefWindow/2
        % Dr Abdul Gaffar: Reference window is not full of data, so cannot reliabilty perform detections for these CUTs  
        gCA = nan;  % Dr Abdul Gaffar
        
    elseif CUT > RefWindow/2 && CUT < length(signal) - RefWindow/2
       LaggingWindow = sum(DataAfterPowerLawDetector( (CUT-RefWindow/2):(CUT-guardCells/2))); 
       LeadingWindow = sum(DataAfterPowerLawDetector( (CUT+guardCells/2):(CUT+RefWindow/2))); 
       gCA = (LaggingWindow + LeadingWindow);
       
    elseif CUT >= length(signal) - RefWindow/2
       % Dr Abdul Gaffar: Reference window is not full of data, so cannot reliabilty perform detections for these CUTs
        gCA = nan;
       
    else
        print('error')
    end
            
    aCA = PFA^(-1/RefWindow)-1; %this is a scaling factor
    
    T_CA_CFAR(CUT) = aCA*gCA;  %threshold value
end

%GOCA_CFAR
T_GOCA_CFAR  = zeros([1  length(signal)]);  %initialise an array for threshold values

for CUT = 1: length(signal)
    if CUT <= RefWindow/2
        % Dr Abdul Gaffar: Reference window is not full of data, so cannot reliabilty perform detections for these CUTs  
        gGO = nan;  % Dr Abdul Gaffar
        
    elseif CUT > RefWindow/2 && CUT < length(signal) - RefWindow/2
        LaggingWindow = sum(DataAfterPowerLawDetector( (CUT-RefWindow/2):(CUT-guardCells/2))); 
        LeadingWindow = sum(DataAfterPowerLawDetector( (CUT+guardCells/2):(CUT+RefWindow/2))); 
        %gCA = LaggingWindow + LeadingWindow
        gGO = max(LeadingWindow, LaggingWindow);
       
    elseif CUT >= length(signal) - RefWindow/2
       % Dr Abdul Gaffar: Reference window is not full of data, so cannot reliabilty perform detections for these CUTs
        gGO = nan;
       
    else
        print('error')
    end
    
    T_GOCA_CFAR(CUT) = aGO*gGO;  %threshold value
end

%SOCA_CFAR
T_SOCA_CFAR  = zeros([1  length(signal)]);  %initialise an array for threshold values

for CUT = 1: length(signal)
    if CUT <= RefWindow/2
        % Dr Abdul Gaffar: Reference window is not full of data, so cannot reliabilty perform detections for these CUTs  
        gSO = nan;  % Dr Abdul Gaffar
        
    elseif CUT > RefWindow/2 && CUT < length(signal) - RefWindow/2
        LaggingWindow = sum(DataAfterPowerLawDetector( (CUT-RefWindow/2):(CUT-guardCells/2))); 
        LeadingWindow = sum(DataAfterPowerLawDetector( (CUT+guardCells/2):(CUT+RefWindow/2))); 
        %gCA = LaggingWindow + LeadingWindow
        gSO = min(LeadingWindow, LaggingWindow);
       
    elseif CUT >= length(signal) - RefWindow/2
       % Dr Abdul Gaffar: Reference window is not full of data, so cannot reliabilty perform detections for these CUTs
        gSO = nan;
       
    else
        print('error')
    end
            
    
    T_SOCA_CFAR(CUT) = aSO*gSO;  %threshold value
end

%OS_CFAR
pos = round(0.75*RefWindow); %Grabs the position of the statistic in the array    
alpha_values = 1:1:500;
Pfa_values = factorial(RefWindow)*factorial(alpha_values+RefWindow-pos)./(factorial(RefWindow-pos)*factorial(alpha_values+RefWindow));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*PFA));
alpha = alpha_values(ind);

T_OS_CFAR  = zeros([1  length(signal)]);  %initialise an array for threshold values

for CUT = 1: length(signal)
    if CUT <= RefWindow/2
        gOS = nan;
        
    elseif CUT > RefWindow/2 && CUT < length(signal) - RefWindow/2
        LaggingWindow = DataAfterPowerLawDetector( (CUT-RefWindow/2):(CUT-guardCells/2)); 
        LeadingWindow = DataAfterPowerLawDetector( (CUT+guardCells/2):(CUT+RefWindow/2)); 
        window = [LeadingWindow, LaggingWindow];
        sortedWindow = sort(window);
        gOS = sortedWindow(pos);
    
    elseif CUT >= length(signal) - RefWindow/2
        gOS = nan; 
       
    else
        print('error')
    end
    
    
    T_OS_CFAR(CUT) = alpha*gOS;  %threshold value
end

for i = 1:1:length(DataAfterPowerLawDetector)
    DataAfterPowerLawDetector(i) = 10*log10(DataAfterPowerLawDetector(i)/10);
end
for i = 1:1:length(T_CA_CFAR)
    T_CA_CFAR(i) = 10*log10(T_CA_CFAR(i)/10);
end
for i = 1:1:length(T_GOCA_CFAR)
    T_GOCA_CFAR(i) = 10*log10(T_GOCA_CFAR(i)/10);
end
for i = 1:1:length(T_SOCA_CFAR)
    T_SOCA_CFAR(i) = 10*log10(T_SOCA_CFAR(i)/10);
end
for i = 1:1:length(T_OS_CFAR)
    T_OS_CFAR(i) = 10*log10(T_OS_CFAR(i)/10);
end
% Plot subset of data

Subset = 500;                                                  % Dr Abdul Gaffar 
%DataAfterPowerLawDetector(201) = 5;

fig4 = figure(4);
ax4 = axes('Parent', fig4);
plot(ax4, t(1:Subset), DataAfterPowerLawDetector(1:Subset))
title('Threshold vs Range')
hold on
plot(ax4, t(1:Subset), T_CA_CFAR(1:Subset))
%plot(ax4, t(1:Subset), T_GOCA_CFAR(1:Subset))
%plot(ax4, t(1:Subset), T_SOCA_CFAR(1:Subset))
%plot(ax4, t(1:Subset), T_OS_CFAR(1:Subset))
legend('Signal (Clutter and Noise)' ,'CA-CFAR', 'GOCA-CFAR', 'SOCA-CFAR', 'OS-CFAR');
xlabel('Range Bins');
ylabel('Threshold value (dB)');
hold off