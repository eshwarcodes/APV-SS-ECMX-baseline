function [panel_reflected_front, panel_reflected_rear] = computePanelReflectedSpectral( ...
    numPanels, groundX,Front_Irradiance,Rear_Irradiance,viewFactor_FrontPVGround,viewFactor_RearPVGround,wavelengthDir,Rloss_Albedo_front_interp,...
   Rloss_Albedo_back_interp)

[Ny, Nx] = size(groundX);
nLambda = length(wavelengthDir);

% Initialize output as a 200x200 cell array, each cell with a 1882x1 vector
panel_reflected_front = cell(Ny, Nx);
panel_reflected_rear = cell(Ny, Nx);
for iy = 1:Ny
    for ix = 1:Nx
        panel_reflected_front{iy, ix} = zeros(nLambda, 1);
        panel_reflected_rear{iy, ix} = zeros(nLambda, 1);
    end
end

% Loop through each panel
%VF_sum_front = zeros(Ny, Nx);
%VF_sum_rear = zeros(Ny, Nx);

for i = 1:numPanels
    % 1882×1 reflected spectra per side
    R_front = Front_Irradiance{i}.*Rloss_Albedo_front_interp;

    R_back = Rear_Irradiance{i}.*Rloss_Albedo_back_interp;
        % Irradiance [W/m²] at ground point
        % Eground = VF_panels{1} * Epanel;
        % 
        % % Total power
        % Eground_total = sum(Eground(:)) * A_cell;
        % disp(['Ground total = ', num2str(Eground_total), ' W']);
        % Eground  = Eground.*A_cell;

    % View factor for this panel: 200x200
    VF1 = viewFactor_FrontPVGround{i};
    VF2 = viewFactor_RearPVGround{i};

    % Add reflected spectra to each ground cell
    for iy = 1:Ny
        for ix = 1:Nx
            if VF1(iy,ix) > 0
                panel_reflected_front{iy, ix} = panel_reflected_front{iy, ix} + (VF1(iy,ix) *R_front);
                %VF_sum_front(iy, ix) = VF_sum_front(iy, ix) + VF1(iy,ix);
            end
            if VF2(iy,ix) > 0
                panel_reflected_rear{iy, ix} = panel_reflected_rear{iy, ix} + (VF2(iy,ix)*R_back);
                %VF_sum_rear(iy, ix) = VF_sum_rear(iy, ix) + VF2(iy,ix);
            end
        end
    end
end
% for iy = 1:Ny
%     for ix = 1:Nx
%         if VF_sum_front(iy, ix) ~= 0
%             panel_reflected_front{iy, ix} = panel_reflected_front{iy, ix} ./ VF_sum_front(iy, ix);
%         else
%             % If VF_sum_front is zero, set the cell to a zeros vector
%             panel_reflected_front{iy, ix} = zeros(nLambda, 1);
%         end
% 
%         if VF_sum_rear(iy, ix) ~= 0
%             panel_reflected_rear{iy, ix} = panel_reflected_rear{iy, ix} ./ VF_sum_rear(iy, ix);
%         else
%             panel_reflected_rear{iy, ix} = zeros(nLambda, 1);
%         end
%     end
% end
end
