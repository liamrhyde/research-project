%Controller

%parameters
PFA = 10^-3;    %Probability of false alarm
referenceCells = 16; %Size of window in question
guardCells = 1; %number of guard cells (on either side, so 1 = 1 cell each side)
FFTwindow = 512;

PFAdoppler = 10^-4;
referenceCellsdoppler = 24;
guardCellsdoppler = 1;
% 
% PFA = 10^-3;    %Probability of false alarm
% referenceCells = 32; %Size of window in question
% guardCells = 1; %number of guard cells (on either side, so 1 = 1 cell each side)
% FFTwindow = 512;
% 
% PFAdoppler = 10^-3;
% referenceCellsdoppler = 32;
% guardCellsdoppler = 1;

sizeOfData = size(RangeProfiles_AfterEqNotch);
xSize = sizeOfData(2);
ySize = sizeOfData(1);
k = mod(ySize, FFTwindow);
newSize = ySize - k;
numberOfIterations = newSize/FFTwindow-1; %IMPORTANT, you lose some end data as data cant be used for difference reference windows (6316 is not a factor of 512, so lose data in 172 rows)

trueDetectionArray = zeros(ySize, xSize);
combined_RD_SingleTest_DetectionArray = [];
CA_CFAR_DetectionArray = CA_CFAR_Plottingdetections_1D(RangeProfiles_AfterEqNotch, PFA, referenceCells, guardCells);
RDcolumnsWithDetections = [];
RTcolumnsWithDetections = [];

windowStart = 1;
windowEnd = FFTwindow;
for iter = 1:1:numberOfIterations

    RD_SingleTest_DetectionArray = Range_Doppler_Plotter(windowStart,windowEnd,RangeProfiles_AfterEqNotch,PRF_Hz, PFAdoppler, referenceCellsdoppler, guardCellsdoppler); %can only accept a max 512 gap
    for i = 1:1:FFTwindow 
        for j = 1:1:xSize
            if RD_SingleTest_DetectionArray(i,j) > 0 
                RDcolumnsWithDetections = [RDcolumnsWithDetections, j];
            end    
        end    
    end


    windowStart = 1+windowEnd;
    windowEnd = windowEnd+FFTwindow;
end
RDcolumnsWithDetections2 = unique(RDcolumnsWithDetections);

for i = 1:1:ySize
    for j = 1:1:xSize
        if CA_CFAR_DetectionArray(i,j)>0
            Lia = ismember(j, RDcolumnsWithDetections2);
            if Lia == 1
                trueDetectionArray(:,j) = CA_CFAR_DetectionArray(:,j);
            end
        end
    end
end

% Plot Range Profiles
fontsize1 = 12;
clims = [-40 0];

% Normalise data to have a peak of 0dB or 1 in linear scale
[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));

% Plot range lines
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);
ylabel('Number of pulses','fontsize',fontsize1);
title('Range lines: CA-CFAR Dual Detector','fontsize',fontsize1);
hold on

for i = 1:1:ySize
    for j = 1:1:xSize
        if trueDetectionArray(i,j) > 0
            text(j,i,'X');
        end
    end
end
hold off
