% -------- Helper function for clamping --------
function out = clamp01(x)
    out = min(max(x, 0), 1);
end