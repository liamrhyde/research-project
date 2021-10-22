function Range_Doppler = Range_Doppler_Plotter(startRow, endRow, RangeProfiles_AfterEqNotch,PRF_Hz, PFA, referenceCells, guardCells)
     %Range Doppler Mapper does all the processing to convert to range
     %doppler map, it can draw a map and map the detections, only accpets
     %fft windows of up tp 512 rows. 
     
     %This passes the data to the dector, column by column (each colum is
     %converted to a row for the detector) and recieves a column back from
     %the detector to create a 2d array of 0s and detections. It maps the
     %data to the colour image (note, in the mapping it converts to doppler
     %velocity (i dunno, it just works

    % Variables (Dr Abdul Gaffar)

    c = 299792458;               
    fc = 10e9;                   %  (Dr Abdul Gaffar) Radar's centre frequency  

    % Processing
    RangeProfiles_AfterEqNotch_original = RangeProfiles_AfterEqNotch;
    RangeProfiles_AfterEqNotch_NewData = RangeProfiles_AfterEqNotch_original(startRow:endRow,:); % subset of data to obtain the Range-Doppler map

    % Plot range lines
    fontsize1 = 12;
    clims = [-40 0];


    %Getting Values
    sizeOfData = size(RangeProfiles_AfterEqNotch_NewData);
    ySize = sizeOfData(1);
    xSize = sizeOfData(2);





    %RangeProfiles_AfterEqNotch(:,1);
    Doppler_HRR_profiles = fftshift(fft(RangeProfiles_AfterEqNotch_NewData,[],1), 1);

    %Detector (CA_CFAR version)
    sizeOfData = size(Doppler_HRR_profiles);
    noColumn = sizeOfData(2);
    noRow = sizeOfData (1);
    detectionArray = [];
    
    for j = 1:1:noColumn
        tranpsoseDoppler_HRR_profiles = Doppler_HRR_profiles(:,j).'; %Tranposes each column (there are 54 columns of 512 rows) to a row of 512 columns
        passingArray = tranpsoseDoppler_HRR_profiles;
        detectionArray = [detectionArray, Range_Doppler_Detector_CA_CFAR(PFA, referenceCells, guardCells, noRow, passingArray)]; 
    end
    
    %detectionArrayTrans = detectionArray.';
    Range_Doppler = detectionArray;
%     % Plot Range Profiles
%     fontsize1 = 12;
%     clims = [-40 0];
% 
%     % Normalise data to have a peak of 0dB or 1 in linear scale
%     [MaxRangeLine MaxIdx] = max(max(abs(Doppler_HRR_profiles)));
% 
%     %Convert to range(m)
%     sizeOfData2 = size(Doppler_HRR_profiles);
%     noColumn = sizeOfData2(2);
% 
%     xSize = 1:1:noColumn;
% 
%     
%     % Plot range lines
%     DopplerFrequency = (PRF_Hz/ySize)*(-ySize/2:(ySize/2-1));
%     % Wavelength = c/PRF_Hz;                                      
%     Wavelength = c/fc;                                      %  (Dr Abdul Gaffar)
%     DopplerVelocity = (DopplerFrequency * Wavelength)/2;
%     
%     figure; axes('fontsize',fontsize1);
%     imagesc(xSize, DopplerVelocity, 20*log10(abs(Doppler_HRR_profiles)./MaxRangeLine),clims);    % Dr Abdul Gaffar
%     colorbar;
%     xlabel('Range (bins)','fontsize',fontsize1);
% 
%     ylabel('Velocity (m/s)','fontsize',fontsize1);
%     
%     title('Range Doppler Map','fontsize',fontsize1);
%     
%     hold on
%     
%     for i = 1:1:noRow
%         for j = 1:1:noColumn
%             if detectionArray(i,j) > 0;
%                 text(j,DopplerVelocity(i),'X');
%             end
%         end
%     end
%     hold off
end