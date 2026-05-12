function sweep = load_bragg_table(filepath, energy_range)
% LOAD_BRAGG_TABLE  Load a (Energy, alpha) Bragg-peak table and populate
%                   sweep.values and sweep.alpha_deg for use with run_rcwa.
%
% The file is tab- or comma-separated with a header row.
% Required columns (case-insensitive): Energy, alpha
% Any additional columns (Efficiency, beta, Cff etc.) are ignored.
%
% USAGE
%   sweep = load_bragg_table('bragg_table.csv')
%       Uses every row in the file as-is.
%
%   sweep = load_bragg_table('bragg_table.csv', 1000:100:2000)
%       Interpolates alpha onto the supplied energy grid.
%       Energies outside the file's range error immediately.
%
% OUTPUT
%   sweep.type      = 'bragg'
%   sweep.values    = energy vector (eV)
%   sweep.alpha_deg = paired grazing angle vector (deg), same length

sweep = struct();

if exist(filepath, 'file') ~= 2
    error('load_bragg_table: file not found: %s', filepath);
end

% --- detect delimiter --------------------------------------------------------
fid = fopen(filepath, 'r');
hdr = strtrim(fgetl(fid));
fclose(fid);

if ~isempty(strfind(hdr, sprintf('\t')))
    delim = '\t';
else
    delim = ',';
end

% --- find required columns ---------------------------------------------------
cols       = strtrim(strsplit(hdr, sprintf(delim)));
cols_lower = lower(cols);

idx_E = find(strcmp(cols_lower, 'energy'), 1);
idx_a = find(strcmp(cols_lower, 'alpha'),  1);

if isempty(idx_E)
    error('load_bragg_table: no ''Energy'' column found in %s', filepath);
end
if isempty(idx_a)
    error('load_bragg_table: no ''alpha'' column found in %s', filepath);
end

% --- read numeric data -------------------------------------------------------
data    = dlmread(filepath, sprintf(delim), 1, 0);
tbl_E   = data(:, idx_E)';   % row vectors
tbl_a   = data(:, idx_a)';

% --- apply energy range if supplied ------------------------------------------
if nargin >= 2 && ~isempty(energy_range)

    E_min = min(tbl_E);
    E_max = max(tbl_E);

    out_of_range = energy_range(energy_range < E_min | energy_range > E_max);
    if ~isempty(out_of_range)
        error(['load_bragg_table: requested energies [%s] eV are outside ' ...
               'table range %.1f – %.1f eV'], ...
               num2str(out_of_range, '%.1f '), E_min, E_max);
    end

    % interpolate alpha onto the requested energy grid
    alpha_interp = interp1(tbl_E, tbl_a, energy_range, 'linear');

    sweep.values    = energy_range;
    sweep.alpha_deg = alpha_interp;

    fprintf('load_bragg_table: interpolated %d points onto requested grid\n', ...
        numel(sweep.values));
else
    sweep.values    = tbl_E;
    sweep.alpha_deg = tbl_a;

    fprintf('load_bragg_table: loaded %d (energy, alpha) pairs from file\n', ...
        numel(sweep.values));
end

sweep.type = 'bragg';

fprintf('  Energy range : %.1f – %.1f eV\n',  min(sweep.values), max(sweep.values));
fprintf('  Alpha range  : %.3f – %.3f deg\n', min(sweep.alpha_deg), max(sweep.alpha_deg));

end