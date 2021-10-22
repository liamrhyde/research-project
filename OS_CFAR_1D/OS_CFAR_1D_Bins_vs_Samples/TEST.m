% Input variables 

PFA = 10^-3;    %Probability of false alarm
RefWindow = 32; %total window size (is divided in 2 for leading and lagging)
guardCells = 2; %total number of guard cells (is divided in 2 for leading and lagging)(must be even number greater than 0)
referenceCells = RefWindow;
N = referenceCells;

pos = round(0.75*referenceCells); %Grabs the position of the statistic in the array    
alpha_values = 1:1:500;
Pfa_values = factorial(N)*factorial(alpha_values+N-pos)./(factorial(N-pos)*factorial(alpha_values+N));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*PFA));
alpha = alpha_values(ind);

signal = RangeProfiles_AfterEqNotch(380,:);

dataSize = length(signal);
threshold = OSCFAR_Detector_1D(PFA, referenceCells, guardCells, dataSize, signal, alpha, pos);
t = 1:1:length(signal);

DataAfterPowerLawDetector = abs(signal).^2; %realising signal power

TCA  = zeros([1  length(signal)]);  %initialise an array for threshold values

for CUT = 1: length(signal)
    if CUT <= referenceCells/2
        gCA = nan;

    elseif CUT > referenceCells/2 && CUT < length(signal) - referenceCells/2
        LaggingWindow = DataAfterPowerLawDetector( (CUT-referenceCells/2):(CUT-guardCells/2)); 
        LeadingWindow = DataAfterPowerLawDetector( (CUT+guardCells/2):(CUT+referenceCells/2)); 
        window = [LeadingWindow, LaggingWindow];
        sortedWindow = sort(window);
        gCA = sortedWindow(pos);

    elseif CUT >= length(signal) - referenceCells/2
        gCA = nan; 

    else
        print('error')
    end


    TCA(CUT) = alpha*gCA;  %threshold value
end
    

fig4 = figure(4);
ax4 = axes('Parent', fig4);
plot(ax4, t, DataAfterPowerLawDetector)
title('Signal Power and Threshold vs Bin')
hold on
plot(ax4, t, TCA)
legend('Signal' , 'Threshold');
xlabel('Bin');
ylabel('Power');
hold off