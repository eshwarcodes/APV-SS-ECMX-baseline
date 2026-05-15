function [IB_Front_inten,ID_Iso_Front_inten,ID_Cir_Front_inten,... 
        ID_Hor_Front_inten,IB_Back_inten,ID_Iso_Rear_inten,...
        ID_Cir_Rear_inten,ID_Hor_Rear_inten,...
        IB_HM1_Front_inten,ID_HM1_Iso_Front_inten,ID_HM1_Cir_Front_inten,... 
        ID_HM1_Hor_Front_inten,IB_HM1_Back_inten,ID_HM1_Iso_Rear_inten,...
        ID_HM1_Cir_Rear_inten,ID_HM1_Hor_Rear_inten,...
        IB_HM2_Front_inten,ID_HM2_Iso_Front_inten,ID_HM2_Cir_Front_inten,... 
        ID_HM2_Hor_Front_inten,IB_HM2_Back_inten,ID_HM2_Iso_Rear_inten,...
        ID_HM2_Cir_Rear_inten,ID_HM2_Hor_Rear_inten] = ...
computePanelHMIntensity(numPanels,IB_Front,IB_Front_HM1,IB_Front_HM2,...
        ID_Iso_Front,...
        ID_Iso_Front_HM1,ID_Iso_Front_HM2,...
        ID_Cir_Front,ID_Cir_Front_HM1,...
        ID_Cir_Front_HM2,ID_Hor_Front,...
        ID_Hor_Front_HM1,ID_Hor_Front_HM2,...
        IB_Back,IB_Back_HM1,...
        IB_Back_HM2,ID_Iso_Rear,...
        ID_Iso_Rear_HM1,ID_Iso_Rear_HM2,...
        ID_Cir_Rear,ID_Cir_Rear_HM1,...
        ID_Cir_Rear_HM2,...
        ID_Hor_Rear,ID_Hor_Rear_HM1,...
        ID_Hor_Rear_HM2,wavelengthDir) 

IB_Front_inten = zeros(numPanels,1);
IB_HM1_Front_inten = zeros(numPanels,1);
IB_HM2_Front_inten = zeros(numPanels,1);
ID_Iso_Front_inten = zeros(numPanels,1);
ID_HM1_Iso_Front_inten = zeros(numPanels,1);
ID_HM2_Iso_Front_inten = zeros(numPanels,1);
ID_Cir_Front_inten = zeros(numPanels,1);
ID_HM1_Cir_Front_inten = zeros(numPanels,1);
ID_HM2_Cir_Front_inten = zeros(numPanels,1);
ID_Hor_Front_inten = zeros(numPanels,1);
ID_HM1_Hor_Front_inten = zeros(numPanels,1);
ID_HM2_Hor_Front_inten = zeros(numPanels,1);
IB_Back_inten = zeros(numPanels,1);
IB_HM1_Back_inten = zeros(numPanels,1);
IB_HM2_Back_inten = zeros(numPanels,1);
ID_Iso_Rear_inten = zeros(numPanels,1);
ID_HM1_Iso_Rear_inten = zeros(numPanels,1);
ID_HM2_Iso_Rear_inten = zeros(numPanels,1);
ID_Cir_Rear_inten = zeros(numPanels,1);
ID_HM1_Cir_Rear_inten = zeros(numPanels,1);
ID_HM2_Cir_Rear_inten = zeros(numPanels,1);
ID_Hor_Rear_inten = zeros(numPanels,1);
ID_HM1_Hor_Rear_inten = zeros(numPanels,1);
ID_HM2_Hor_Rear_inten = zeros(numPanels,1);

for k=1:numPanels 
        
    IB_Front_inten(k,1) = trapz(wavelengthDir,IB_Front{k});
    IB_HM1_Front_inten(k,1) = trapz(wavelengthDir,IB_Front_HM1{k});
    IB_HM2_Front_inten(k,1) = trapz(wavelengthDir,IB_Front_HM2{k});    
    
    ID_Iso_Front_inten(k,1) = trapz(wavelengthDir,ID_Iso_Front);
    ID_HM1_Iso_Front_inten(k,1)  = trapz(wavelengthDir,ID_Iso_Front_HM1);
    ID_HM2_Iso_Front_inten(k,1)  = trapz(wavelengthDir,ID_Iso_Front_HM2);
    
    ID_Cir_Front_inten(k,1) = trapz(wavelengthDir,ID_Cir_Front);
    ID_HM1_Cir_Front_inten(k,1)  = trapz(wavelengthDir,ID_Cir_Front_HM1);
    ID_HM2_Cir_Front_inten(k,1)  = trapz(wavelengthDir,ID_Cir_Front_HM2);  

    ID_Hor_Front_inten(k,1) = trapz(wavelengthDir,ID_Hor_Front);
    ID_HM1_Hor_Front_inten(k,1)  = trapz(wavelengthDir,ID_Hor_Front_HM1);
    ID_HM2_Hor_Front_inten(k,1)  = trapz(wavelengthDir,ID_Hor_Front_HM2);
        
    IB_Back_inten(k,1) = trapz(wavelengthDir,IB_Back{k});
    IB_HM1_Back_inten(k,1) = trapz(wavelengthDir,IB_Back_HM1{k});
    IB_HM2_Back_inten(k,1) = trapz(wavelengthDir,IB_Back_HM2{k}); 

    ID_Iso_Rear_inten(k,1) = trapz(wavelengthDir,ID_Iso_Rear);
    ID_HM1_Iso_Rear_inten(k,1)  = trapz(wavelengthDir,ID_Iso_Rear_HM1);
    ID_HM2_Iso_Rear_inten(k,1)  = trapz(wavelengthDir,ID_Iso_Rear_HM2);
        
    ID_Cir_Rear_inten(k,1) = trapz(wavelengthDir,ID_Cir_Rear);
    ID_HM1_Cir_Rear_inten(k,1)  = trapz(wavelengthDir,ID_Cir_Rear_HM1);
    ID_HM2_Cir_Rear_inten(k,1)  = trapz(wavelengthDir,ID_Cir_Rear_HM2);
        
    ID_Hor_Rear_inten(k,1) = trapz(wavelengthDir,ID_Hor_Rear);
    ID_HM1_Hor_Rear_inten(k,1)  = trapz(wavelengthDir,ID_Hor_Rear_HM1);
    ID_HM2_Hor_Rear_inten(k,1)  = trapz(wavelengthDir,ID_Hor_Rear_HM2);
end