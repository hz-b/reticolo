%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     HELPER FUNCTION: Load Custom Refractive Index File (CXRO Format)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nData, filename] = load_custom_refractive_index(filename, reference_dir)

    
    % Construct full path
    if isempty(reference_dir)
        filepath = filename;
    else
        filepath = fullfile(reference_dir, filename);
    end
    
    % Check if file exists
    if ~exist(filepath, 'file')
        error('Refractive index file not found: %s', filepath);
    end
    
    filename = filepath;
    
    % Read the file (skip comment lines starting with #)
    data = [];
    fid = fopen(filepath, 'r');
    if fid == -1
        error('Cannot open file: %s', filepath);
    end
    
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line) && line(1) ~= '#'
            % Parse space/tab separated values
            parts = strsplit(strtrim(line));
            if length(parts) >= 3
                try
                    energy = str2double(parts{1});
                    delta = str2double(parts{2});
                    beta = str2double(parts{3});
                    if ~isnan(energy) && ~isnan(delta) && ~isnan(beta)
                        data = [data; energy, delta, beta];
                    end
                catch
                    % Skip malformed lines
                    continue;
                end
            end
        end
    end
    fclose(fid);
    
    if isempty(data)
        error('No valid data found in file: %s', filepath);
    end
    
    nData = data;
    disp(['Loaded ', num2str(size(nData, 1)), ' data points from ', filename]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     END OF FILE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%