function VF_Panel_to_HM = computeReciprocalVF_PanelToHM( ...
    VF_HM_to_Panel, A_HM,A_PV)

% VF_HM_to_Panel: {panel, face} → {1x2 cell} [HM1, HM2]
% Output:
% VF_Panel_to_HM: same structure, reversed direction via reciprocity

numPanels = size(VF_HM_to_Panel, 1);


% Initialize output
VF_Panel_to_HM = cell(numPanels, 2);  % [panel, face]

for idx = 1:numPanels
    for face = 1:2  % 1 = front, 2 = back
        VF_Panel_to_HM{idx, face} = cell(1, 2);  % [HM1, HM2]
        for h = 1:2
            VF_Panel_to_HM{idx, face}{h} = ...
                (A_HM / A_PV) * VF_HM_to_Panel{idx, face}{h};
        end
    end
end
end
