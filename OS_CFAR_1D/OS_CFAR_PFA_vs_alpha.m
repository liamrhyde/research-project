% Input variables 
Length = 500000; %size of the data (how many samples of noise)

RefWindow = 32; %total window size (is divided in 2 for leading and lagging)

guardCells = 2; %total number of guard cells (is divided in 2 for leading and lagging)(must be even number greater than 0)
alpha_range = 4:1:20;

pos = round(0.75*RefWindow); %Grabs the position of the statistic in the array

Pfa_array = [10^-2, 5*10^-3, 10^-3, 5*10^-4, 10^-4, 5*10^-5, 10^-5,  5*10^-6, 10^-6,];
Pfa_average = [];
Pfa = [];


signal = zeros([1  Length]); %signal initialisation
t = 1:1:length(signal); %time of simulation
signal = normrnd(0,1,1,Length) + 1i*normrnd(0,1,1,Length);  %complex white noise generation

DataAfterPowerLawDetector = abs(signal).^2; %realising signal power

TCA  = zeros([1  length(signal)]);  %initialise an array for threshold values



for alpha  = alpha_range
    
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
    Pfa = [Pfa;PFA_simulation];
   
    Pfa_average = [Pfa_average; factorial(RefWindow)*factorial(alpha+RefWindow-pos)/(factorial(RefWindow-pos)*factorial(alpha+RefWindow))];
   
end;
close all
figure;
semilogy(alpha_range,Pfa);
hold on;
semilogy(alpha_range,Pfa_average);
grid on;
xlabel('Alpha Value');
ylabel('Pfa');
title('OS-CFAR Pfa vs alpha')
legend('Simulated Pfa','Expected Pfa')
