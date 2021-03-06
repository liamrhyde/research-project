%IMPORTANT
%This program plots range (bins) vs Samples, placing an X over each
%detection, in 1D (each row is analysed)


%parameters
PFA = 10^-5;    %Probability of false alarm NOTE: it has a lower PFA than CA_CFAR and GOCA_CFAR
referenceCells = 16; %total window size (is divided in 2 for leading and lagging)
guardCells = 2; %total number of guard cells (is divided in 2 for leading and lagging)(must be even number greater than 0)
aCA = SOCA_CFAR_Alpha(PFA, referenceCells);

%Getting sizes of data
sizeOfData = size(RangeProfiles_AfterEqNotch);
noColumn = sizeOfData(2);
dataSize = noColumn;
noRow = sizeOfData (1);
detectionArray = [];

    

%Gets 1 row of data from the RangeProfile and sends it to CACFAR, which
%then returns a single row containing only the values that were detected as
%targets and appends it to the detection array. Detection array is then
%sent to data extraction to mark each detection with an x (test with PFA =
%10^-1, refcells - 9, guard cells = 0 and rangeprofile columns 1 - 24)


for i = 1:1:noRow;    
    passingArray = RangeProfiles_AfterEqNotch(i,:);
    
    detectionArray = [detectionArray; SOCACFAR_Detector_1D(PFA, referenceCells, guardCells, dataSize, passingArray, aCA)];
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
title('Range lines: after Eq Notch','fontsize',fontsize1);
hold on

for i = 1:1:noRow;
    for j = 1:1:noColumn;
        if detectionArray(i,j) > 0;
            text(j,i,'X');
        end
    end
end
hold off
