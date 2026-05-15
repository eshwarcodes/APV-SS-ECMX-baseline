function status = stopAfterMinutes(t, y, flag)
    persistent startTime maxDuration

    if isempty(flag)
        elapsed = toc(startTime);
        if elapsed > maxDuration
            status = 1; % Stop the integration
        else
            status = 0; % Continue
        end
    elseif strcmp(flag, 'init')
        startTime = tic;
        maxDuration = 60 * 1; % ← change this to your desired duration in seconds (e.g., 10 min)
        status = 0;
    elseif strcmp(flag, 'done')
        clear startTime maxDuration
        status = 0;
    end
end
