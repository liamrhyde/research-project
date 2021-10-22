%IMPORTANT
%This program plots range (bins) vs Samples, placing an X over each
%detection, in 1D (each row is analysed)


%parameters
% PFA = 10^-4;    %Probability of false alarm NOTE: it has a lower PFA than CA_CFAR
% referenceCells = 18; %total window size (is divided in 2 for leading and lagging)
% guardCells = 2; %total number of guard cells (is divided in 2 for leading and lagging)(must be even number greater than 0)


PFA = 10^-3;    %Probability of false alarm NOTE: it has a lower PFA than CA_CFAR
referenceCells = 32; %total window size (is divided in 2 for leading and lagging)
guardCells = 2; %total number of guard cells (is divided in 2 for leading and lagging)(must be even number greater than 0)

aCA = GOCA_CFAR_Alpha(PFA, referenceCells);

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


for i = 1:1:noRow    
    passingArray = RangeProfiles_AfterEqNotch(i,:);
    
    detectionArray = [detectionArray; GOCACFAR_Detector_1D(PFA, referenceCells, guardCells, dataSize, passingArray, aCA)];
end

% Plot Range Profiles
fontsize1 = 12;
clims = [-40 0];

% Normalise data to have a peak of 0dB or 1 in linear scale
[MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));

%Convert to range(m)
sizeOfData = size(RangeProfiles_AfterEqNotch);
noColumn = sizeOfData(2);

%need to get number of columns to nearest 5 or 10
k = mod(noColumn, 5);
if k > 0;
    newNoColumn = noColumn - k;
else
    newNoColumn = noColumn;
end
if newNoColumn == 800;
    xStep = newNoColumn/8;
elseif newNoColumn >= 100;
    k2 = mod(newNoColumn,100);
    newNoColumn = newNoColumn-k2;
    xStep = round(newNoColumn/10);
else
    xStep = round(newNoColumn/10);
end
xSize = 1:1:noColumn;
distanceInMeters = xSize(xStep:xStep:end);
%distanceInMeters = xSize*299792458/(2*Bandwidth_Hz);
distanceInMetersScaled = distanceInMeters*(299792458/(2*Bandwidth_Hz));

%Convert to time(s)
noRow = sizeOfData(1);
PRI = 1/PRF_Hz;

k = mod(noRow, 500);
if k > 0;
    newNoRow = noRow - k;
else
    newNoRow = noRow;
end
if newNoRow == 4000;
    yStep = newNoRow/8;
elseif newNoRow == 6000;
    yStep = newNoRow/6;
elseif newNoRow == 3000;
    yStep = newNoRow/6;
else
    yStep = round(newNoColumn/10);
end
ySize = 1:1:noRow;
timeInSeconds = ySize(yStep:yStep:end);
timeInSecondsScaled = timeInSeconds*PRI;


% Plot range lines
figure; axes('fontsize',fontsize1);
imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
colorbar;
xlabel('Range (meters)','fontsize',fontsize1);
xticklabels(cellstr(num2str(distanceInMetersScaled.')));
ylabel('Time (seconds)','fontsize',fontsize1);
yticklabels(cellstr(num2str(timeInSecondsScaled.')));
title('Range vs Time: GOCA','fontsize',fontsize1);

for i = 1:1:noRow;
    for j = 1:1:noColumn;
        if detectionArray(i,j) > 0;
            text(j,i,'X');
        end
    end
end
hold off
