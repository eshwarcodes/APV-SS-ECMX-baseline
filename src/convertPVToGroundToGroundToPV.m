function [viewFactorCell_GroundToPV,viewFactorCell_GroundToPV_sparse] = convertPVToGroundToGroundToPV( ...
    viewFactorCell_PVToGround, panelWidth, panelHeight, ...
    groundXmin, groundXmax, groundYmin, groundYmax, numGroundX, numGroundY)

% Number of panels
numPanels = numel(viewFactorCell_PVToGround);
viewFactorCell_GroundToPV = cell(numPanels, 1);
viewFactorCell_GroundToPV_sparse = cell(numPanels, 1);
% Area of each panel
A_panel = panelWidth * panelHeight;

% Area of one ground cell
dx = (groundXmax - groundXmin) / numGroundX;
dy = (groundYmax - groundYmin) / numGroundY;
A_ground_cell = dx * dy;

% === Step 1: Apply reciprocity to each panel ===
for i = 1:numPanels
    VF_pv_to_ground = full(viewFactorCell_PVToGround{i});  % 200x200
    viewFactorCell_GroundToPV{i} = VF_pv_to_ground * (A_panel / A_ground_cell);
end
% imagesc(viewFactorCell_GroundToPV{1});
% colorbar;
% title('Unnormalized VF from Ground to Panel 1');
% === Step 2: Normalize per ground cell (so total VF ≤ 1 from each ground point) ===
VF_sum = zeros(numGroundY, numGroundX);  % accumulate total VF from each ground cell

for i = 1:numPanels
    VF_sum = VF_sum + viewFactorCell_GroundToPV{i};
end

% Avoid division by zero using a small epsilon
VF_sum_safe = max(VF_sum, 1e-12);

for i = 1:numPanels
    viewFactorCell_GroundToPV{i} = viewFactorCell_GroundToPV{i} ./ VF_sum_safe;
    viewFactorCell_GroundToPV_sparse{i} = single(viewFactorCell_GroundToPV{i});
end

end
