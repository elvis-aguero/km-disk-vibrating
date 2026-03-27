function ramBytes = getAvailableRAM()
    % GETAVAILABLERAM Returns the available physical RAM in bytes.
    % Works on macOS, Linux, and Windows.
    
    if ispc
        % Windows: Use wmic
        [~, out] = system('wmic OS get FreePhysicalMemory /Value');
        ramBytes = str2double(regexp(out, '\d+', 'match')) * 1024;
    elseif ismac
        % macOS: Use sysctl for total physical memory (simpler proxy)
        % Note: Getting truly "available" memory on macOS is complex via CLI.
        % We use total physical memory as the limit for the gatekeeper.
        [~, out] = system('sysctl hw.memsize');
        ramBytes = str2double(regexp(out, '\d+', 'match'));
    else 
        % Linux: Read MemAvailable from /proc/meminfo
        [status, out] = system('grep MemAvailable /proc/meminfo');
        if status == 0
            ramBytes = str2double(regexp(out, '\d+', 'match')) * 1024;
        else
            % Fallback for older kernels
            [~, out] = system('grep MemFree /proc/meminfo');
            ramBytes = str2double(regexp(out, '\d+', 'match')) * 1024;
        end
    end
    
    % Final fallback if all else fails
    if isnan(ramBytes) || ramBytes == 0
        ramBytes = 8e9; % Default to 8GB
    end
end
