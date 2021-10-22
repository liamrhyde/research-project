function CA = Range_Doppler_Detector_CA_CFAR(PFA, referenceCells, guardCells, noRow, passingArray)
    
    %passingArray is a single column of the columns in the range doppler
    %map, but have been converted to a row for the detector. Each column is being passed 1 at a time, so detector operates one
    %at a time
    
    %returns a detection row, that is converted to a column before sending
    
    DataAfterPowerLawDetector = abs(passingArray).^2; %realising signal power
    TCA  = zeros([1  length(passingArray)]);    %a 1 by noRow array
    
    
    for CUT = 1: length(passingArray)
        if CUT <= referenceCells/2
            gCA = nan;
        
        elseif CUT > referenceCells/2 && CUT < length(passingArray) - referenceCells/2
            LaggingWindow = sum(DataAfterPowerLawDetector( (CUT-referenceCells/2):(CUT-guardCells))); 
            LeadingWindow = sum(DataAfterPowerLawDetector( (CUT+guardCells):(CUT+referenceCells/2))); 
            gCA = LaggingWindow + LeadingWindow;
       
        elseif CUT >= length(passingArray) - referenceCells/2
            gCA = nan; 
       
        else
            print('error')
        end
        
        aCA = PFA^(-1/referenceCells)-1; %this is a scaling factor
        TCA(CUT) = aCA*gCA;  %threshold value
    end
    
    %checking detections
    detections = zeros([1 length(passingArray)]);
    % PFA_simulation = 0;
    NumberOfFalseAlarms = 0; % Dr Abdul Gaffar

    for i = 1:length(passingArray)
        if DataAfterPowerLawDetector(i) >= TCA(i)
            detections(i) = DataAfterPowerLawDetector(i);
            NumberOfFalseAlarms = NumberOfFalseAlarms + 1;
        end
    end
    
    PFA_simulation = NumberOfFalseAlarms/(noRow - referenceCells); 
    PFA_error = abs(((PFA - PFA_simulation)/PFA)*100) ;  
    
                       
    CA = detections.'; %converts back to a column
    
end