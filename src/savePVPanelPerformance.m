function savePVPanelPerformance(filename,panelResults,Front_Irradiance_abs,Rear_Irradiance_abs)
       
        % Save the variable in the specified subfolder
        save(filename, 'panelResults','Front_Irradiance_abs','Rear_Irradiance_abs');
end