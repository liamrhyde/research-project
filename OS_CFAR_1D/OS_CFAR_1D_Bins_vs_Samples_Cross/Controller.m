%Controller

%parameters RANGE - TIME
PFA_T = 10^-3;    %Probability of false alarm
referenceCells_T = 32; %Size of window in question
guardCells_T = 2; %number of guard cells (divided by 2, so 2 = 1 cell each side)

pos_T = round(0.75*referenceCells_T); %Grabs the position of the statistic in the array 

alpha_values = 1:1:500;
Pfa_values = factorial(referenceCells_T)*factorial(alpha_values+referenceCells_T-pos_T)./(factorial(referenceCells_T-pos_T)*factorial(alpha_values+referenceCells_T));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*PFA_T));
alpha_T = alpha_values(ind);

%parameters RANGE - DOPPLER
PFA_D = 10^-3;
referenceCells_D = 32;
guardCells_D = 2;

pos_D = round(0.75*referenceCells_D); %Grabs the position of the statistic in the array   

alpha_values = 1:1:500;
Pfa_values = factorial(referenceCells_D)*factorial(alpha_values+referenceCells_D-pos_D)./(factorial(referenceCells_D-pos_D)*factorial(alpha_values+referenceCells_D));
[val,ind] =  min(abs(Pfa_values - ones(1,length(alpha_values))*PFA_D));
alpha_D = alpha_values(ind);

%Data Size
sizeOfData = size(RangeProfiles_AfterEqNotch);
xSize = sizeOfData(2);
ySize = sizeOfData(1);

%RANGE - TIME DETECTIONS
%OS_CFAR_DetectionArray = OS_CFAR_Plottingdetections_1D(RangeProfiles_AfterEqNotch, PFA_T, referenceCells_T, guardCells_T, alpha_T, pos_T);
%OS_CFAR_DetectionArray is an array of 0's and detections, in the size of
%the original matrix, RangeProfiles_AfterEqNotch.


%CALCULATING WINDOW SIZES TO PASS FOR FFT
FFTwindow = 512;
k = mod(ySize, FFTwindow);
newSize = ySize - k;
numberOfIterations = newSize/FFTwindow-1; %IMPORTANT, you lose some end data as data cant be used for difference reference windows (6316 is not a factor of 512, so lose data in 172 rows)


%PASSING DATA FOR FFT

doppler_to_time_combined_array = [];

windowStart = 1;
windowEnd = FFTwindow;

for iter = 1:1:1

    RD_SingleTest_DetectionArray = Range_Doppler_Plotter(windowStart,windowEnd,RangeProfiles_AfterEqNotch,PRF_Hz, PFA_D, referenceCells_D, guardCells_D, alpha_D, pos_D); %can only accept a max 512 gap
    %doppler_to_time_combined_array = [doppler_to_time_combined_array;RD_SingleTest_DetectionArray];


    windowStart = 1+windowEnd;
    windowEnd = windowEnd+FFTwindow;
end