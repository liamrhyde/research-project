
%% Clear workspace and all figures
clear all;clc;close all

%% Load Data
[FileName,Path2RadarData,filter_index]=uigetfile('','Select Radar Dataset');
RadarData = load([Path2RadarData filesep FileName])

%% Extract Range Profiles before, after Equalisation and after Notch filtering

RangeProfiles_BeforeEq = RadarData.RangeLines_BeforeEq;
RangeProfiles_AfterEq = RadarData.RangeLines_AfterEq;
RangeProfiles_AfterEqNotch = RadarData.RangeLines_AfterEQ_Notch;
%% Extract other radar parameters

 PRF_Hz = RadarData.Info.PRF_Hz;
 Bandwidth_Hz = RadarData.Info.Bandwidth_Hz;
 RangeStart_m = RadarData.Info.RangeStart_m;
 BlindRange_m = RadarData.Info.BlindRange_m;
 
[NumOfPulses,NumOfRangeBins]=size(RangeProfiles_AfterEqNotch);

%% Plot Range Profiles

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
title('Range lines: after Eq Notch','fontsize',fontsize1);
