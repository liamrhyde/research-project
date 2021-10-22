%What the verification does is generate a number of samples of white noise
%with no taget present. This means that 100% of detections are false
%alarms. By setting the false alarm rate and comparing the expect number of
%false alarms with the actual number of false alarms the error can be
%calculated.. error should be less than 10% for it to be correct

close all;
clear all;

% Input variables 
Length = 500000; %size of the data (how many samples of noise)
PFA = 10^-4;    %Probability of false alarm
RefWindow = 16; %total window size (is divided in 2 for leading and lagging)
guardCells = 2; %total number of guard cells (is divided in 2 for leading and lagging)(must be even number greater than 0)

pos = round(0.75*RefWindow); %Grabs the position of the statistic in the array    
alpha_values = 1:1:500;
Pfa_values = factorial(RefWindow)*factorial(alpha_values+RefWindow-pos)./(factorial(RefWindow-pos)*factorial(alpha_values+RefWindow));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*PFA));
alpha = alpha_values(ind)


signal = zeros([1  Length]); %signal initialisation
t = 1:1:length(signal); %time of simulation
signal = normrnd(0,1,1,Length) + 1i*normrnd(0,1,1,Length);  %complex white noise generation

DataAfterPowerLawDetector = abs(signal).^2; %realising signal power

TCA  = zeros([1  length(signal)]);  %initialise an array for threshold values

for CUT = 1: length(signal)
    if CUT <= RefWindow/2
        gCA = nan;
        
    elseif CUT > RefWindow/2 && CUT < length(signal) - RefWindow/2
        LaggingWindow = DataAfterPowerLawDetector( (CUT-RefWindow/2):(CUT-guardCells/2)); 
        LeadingWindow = DataAfterPowerLawDetector( (CUT+guardCells/2):(CUT+RefWindow/2)); 
        window = [LeadingWindow, LaggingWindow];
        sortedWindow = sort(window);
        gCA = sortedWindow(pos);
    
    elseif CUT >= length(signal) - RefWindow/2
        gCA = nan; 
       
    else
        print('error')
    end
    
    
    TCA(CUT) = alpha*gCA;  %threshold value
end
    
%checking detections
detections = zeros([1 length(signal)]);

NumberOfFalseAlarms = 0; % Dr Abdul Gaffar

for i = 1:length(signal)
    if DataAfterPowerLawDetector(i) >= TCA(i)
        detections(i) = DataAfterPowerLawDetector(i);
        NumberOfFalseAlarms = NumberOfFalseAlarms + 1; 
    end
end

PFA_simulation = NumberOfFalseAlarms/(Length - RefWindow); 

%checking false alarms


PFA_error = abs(((PFA - PFA_simulation)/PFA)*100)               % Dr Abdul Gaffar 


% Plot subset of data

Subset = 5000;                                                  % Dr Abdul Gaffar 

fig4 = figure(4);
ax4 = axes('Parent', fig4);
plot(ax4, t(1:Subset), DataAfterPowerLawDetector(1:Subset))
title('Signal Power and Threshold vs Time')
hold on
plot(ax4, t(1:Subset), TCA(1:Subset))
legend('Signal' , 'Threshold');
hold off

fig5 = figure(5);
ax5 = axes('Parent', fig5);
plot(ax5, t(1:Subset), detections(1:Subset))
title('Detections vs Time')