function VF_HM_to_Panel = computeRayBasedHMToPanelVF( ...
    hmCorners, panelCorners, n_panel, n_hm1, n_hm2, ...
    N_rays, maxViewAngle_deg)

numPanels = length(hmCorners);
VF_HM_to_Panel = cell(numPanels, 2);  % (panel idx, face) → each is 1x2 cell [HM1, HM2]

cos_maxView = cosd(maxViewAngle_deg);

for idx = 1:numPanels
    VF_front = cell(1,2);
    VF_back  = cell(1,2);

    for h = 1:2
        if isempty(hmCorners{idx}) || isempty(panelCorners{idx})
            VF_front{h} = 0;
            VF_back{h} = 0;
            continue;
        end

        if h == 1
            n_hm = n_hm1;
        else
            n_hm = n_hm2;
        end

        % --- Hot mirror corners
        C = hmCorners{idx}{h};
        p1 = C(1,:); p2 = C(2,:); p3 = C(3,:); p4 = C(4,:);

        % --- Panel corners
        P = panelCorners{idx};
        n_panel_front = n_panel;
        n_panel_back  = -n_panel;

        % Dimensions
        Lx = norm(p2 - p1);
        Ly = norm(p4 - p1);
        nX = ceil(sqrt(N_rays));
        nY = nX;
        dx = Lx / nX;
        dy = Ly / nY;

        % Axes
        u = (p2 - p1) / norm(p2 - p1);  % width
        v = (p4 - p1) / norm(p4 - p1);  % height

        dA = dx * dy;
        VF_hm_to_front = 0;
        VF_hm_to_back = 0;

        for i = 1:nX
            for j = 1:nY
                sample = p1 + (i - 0.5) * dx * u + (j - 0.5) * dy * v;

                for corner = 1:4
                    Q = P(corner,:);
                    r_vec = Q - sample;
                    r = norm(r_vec);
                    if r < 1e-4, continue; end
                    r_hat = r_vec / r;

                    cos_theta_HM    = dot(n_hm, r_hat);
                    cos_theta_front = dot(n_panel_front, -r_hat);
                    cos_theta_back  = dot(n_panel_back, -r_hat);

                    if cos_theta_HM > 0 && cos_theta_front > 0 && cos_theta_front >= cos_maxView
                        dVF = (cos_theta_HM * cos_theta_front) / (pi * r^2) * dA;
                        VF_hm_to_front = VF_hm_to_front + dVF;
                    end

                    if cos_theta_HM > 0 && cos_theta_back > 0 && cos_theta_back >= cos_maxView
                        dVF = (cos_theta_HM * cos_theta_back) / (pi * r^2) * dA;
                        VF_hm_to_back = VF_hm_to_back + dVF;
                    end
                end
            end
        end

        VF_front{h} = VF_hm_to_front;
        VF_back{h}  = VF_hm_to_back;
    end

    VF_HM_to_Panel{idx,1} = VF_front;  % Front face of panel
    VF_HM_to_Panel{idx,2} = VF_back;   % Rear face of panel
end
end
