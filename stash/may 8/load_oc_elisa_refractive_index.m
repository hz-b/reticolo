function [nData, nfile] = load_oc_elisa_refractive_index(material, reference_dir, density)
%LOAD_OC_ELISA_REFRACTIVE_INDEX Load OC_ELISA refractive index CSV data.
%   The file format is:
%     OC_ELISA_<material>_p<density>.csv
%   with semicolon separators and comma decimal notation.

if nargin < 2 || isempty(reference_dir)
    reference_dir = '';
end
if nargin < 3
    density = [];
end

candidate_dirs = {};
if ~isempty(reference_dir)
    candidate_dirs{end + 1} = reference_dir;
end
candidate_dirs{end + 1} = fullfile(pwd, 'RE_Sample_RF_file');
candidate_dirs{end + 1} = fullfile(pwd, 'RE__Sample_RF_file');

reference_dir = '';
for k = 1:numel(candidate_dirs)
    if exist(candidate_dirs{k}, 'dir') == 7
        reference_dir = candidate_dirs{k};
        break;
    end
end

if isempty(reference_dir)
    error('Reference folder not found. Expected RE_Sample_RF_file or RE__Sample_RF_file.');
end

files = dir(fullfile(reference_dir, sprintf('OC_ELISA_%s_p*.csv', material)));
if isempty(files)
    error('No OC_ELISA refractive index file found for material "%s" in %s.', material, reference_dir);
end

file_names = {files.name};
selected_index = [];

if ~isempty(density)
    density_text = num2str(density, '%.15g');
    exact_name = sprintf('OC_ELISA_%s_p%s.csv', material, density_text);
    exact_index = find(strcmp(file_names, exact_name), 1);
    if ~isempty(exact_index)
        selected_index = exact_index;
    else
        densities = nan(size(file_names));
        for k = 1:numel(file_names)
            tokens = regexp(file_names{k}, '^OC_ELISA_.+_p([0-9]+(?:\.[0-9]+)?)\.csv$', 'tokens', 'once');
            if ~isempty(tokens)
                densities(k) = str2double(tokens{1});
            end
        end
        [~, selected_index] = min(abs(densities - density));
    end
else
    [~, selected_index] = sort(lower(file_names));
    selected_index = selected_index(1);
    if numel(files) > 1
        warning('Multiple refractive index files found for "%s". Using %s. Set a density to disambiguate.', ...
            material, file_names{selected_index});
    end
end

nfile = fullfile(reference_dir, file_names{selected_index});
fid = fopen(nfile, 'r');
if fid < 0
    error('Unable to open refractive index file: %s', nfile);
end
cleanup = onCleanup(@() fclose(fid));

header_line = fgetl(fid); %#ok<NASGU>
raw = textscan(fid, '%s%s%s', 'Delimiter', ';');
if isempty(raw) || isempty(raw{1})
    error('No refractive index data found in file: %s', nfile);
end

energy = str2double(strrep(raw{1}, ',', '.'));
delta = str2double(strrep(raw{2}, ',', '.'));
beta = str2double(strrep(raw{3}, ',', '.'));

nData = [energy, delta, beta];
end
