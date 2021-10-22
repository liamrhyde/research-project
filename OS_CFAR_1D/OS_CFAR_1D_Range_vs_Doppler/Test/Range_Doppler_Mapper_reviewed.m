%Range Doppler Mapper

% Variables (Dr Abdul Gaffar)
 NumOfProfiles2Process = 2048; %  (Dr Abdul Gaffar) to generate one Range-Doppler map
                              %  (Dr Abdul Gaffar) Processing in time (seconds) =  NumOfProfiles2Process*1/PRF_Hz
 c = 299792458;               %  (Dr Abdul Gaffar)  
 fc = 10e9;                   %  (Dr Abdul Gaffar) Radar's centre frequency  
 
 % Processing
 RangeProfiles_AfterEqNotch_original = RangeProfiles_AfterEqNotch;
 RangeProfiles_AfterEqNotch = RangeProfiles_AfterEqNotch_original(1:NumOfProfiles2Process,:); % subset of data to obtain the Range-Doppler map
 
 % Plot range lines
 fontsize1 = 12;
 clims = [-40 0];

 % Normalise data to have a peak of 0dB or 1 in linear scale
 [MaxRangeLine MaxIdx] = max(max(abs(RangeProfiles_AfterEqNotch)));

 figure; axes('fontsize',fontsize1);
 imagesc(20*log10(abs(RangeProfiles_AfterEqNotch)./MaxRangeLine),clims);
 colorbar;
 xlabel('Range (bins)','fontsize',fontsize1);
 ylabel('Time (bins)','fontsize',fontsize1);
 title('Range lines: after Eq Notch','fontsize',fontsize1);
 
%Getting Values
sizeOfData = size(RangeProfiles_AfterEqNotch);
ySize = sizeOfData(1);
xSize = sizeOfData(2);

%parameters
PFA = 10^-2;    %Probability of false alarm
referenceCells = 16; %Size of window in question
guardCells = 0; %number of guard cells (on either side, so 1 = 1 cell each side)

dataAfterFFT = zeros(ySize,xSize);

%RangeProfiles_AfterEqNotch(:,1);
Doppler_HRR_profiles = fftshift(fft(RangeProfiles_AfterEqNotch,[],1), 1);


% Plot Range Profiles
fontsize1 = 12;
clims = [-40 0];

% Normalise data to have a peak of 0dB or 1 in linear scale
[MaxRangeLine MaxIdx] = max(max(abs(Doppler_HRR_profiles)));

%Convert to range(m)
sizeOfData2 = size(Doppler_HRR_profiles);
noColumn = sizeOfData2(2);

%need to get number of columns to nearest 5 or 10
k = mod(noColumn, 5);
if k > 0
    newNoColumn = noColumn - k;
else
    newNoColumn = noColumn;
end
if newNoColumn == 800
    xStep = newNoColumn/8;
elseif newNoColumn >= 100
    k2 = mod(newNoColumn,100);
    newNoColumn = newNoColumn-k2;
    xStep = round(newNoColumn/10);
else
    xStep = round(newNoColumn/10);
end

xSize = 1:1:noColumn;
distanceInMeters = xSize(xStep:xStep:end);
%distanceInMeters = xSize*299792458/(2*Bandwidth_Hz);
distanceInMetersScaled = distanceInMeters*(c/(2*Bandwidth_Hz));


% Plot range lines
DopplerFrequency = (PRF_Hz/ySize)*(-ySize/2:(ySize/2-1))
% Wavelength = c/PRF_Hz;                                      
Wavelength = c/fc;                                      %  (Dr Abdul Gaffar)
DopplerVelocity = (DopplerFrequency * Wavelength)/2;
figure; axes('fontsize',fontsize1);
% imagesc(distanceInMetersScaled,DopplerVelocity,20*log10(abs(Doppler_HRR_profiles)./MaxRangeLine),clims);
imagesc(xSize,DopplerVelocity,20*log10(abs(Doppler_HRR_profiles)./MaxRangeLine),clims);    % Dr Abdul Gaffar
colorbar;
xlabel('Range (bins)','fontsize',fontsize1);

ylabel('Velocity (m/s)','fontsize',fontsize1);

title('Range lines: after Eq Notch','fontsize',fontsize1);