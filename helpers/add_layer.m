function stack = add_layer(stack, material_file, thickness_nm)
% ADD_LAYER  Append one conformal coating layer above the existing stack.
%
% Layers are conformal: each layer follows the shape of the surface below
% it, offset upward by thickness_nm in z everywhere.
%
% INPUTS
%   stack          the stack struct (from build_stack or a previous add_layer)
%   material_file  path to CXRO-format optical constants file
%                  expected columns:  Energy(eV)  Delta  Beta
%                  with header lines starting with non-numeric characters
%   thickness_nm   layer thickness in nm
%
% OUTPUT
%   stack          updated stack with the new layer appended

if thickness_nm <= 0
    error('add_layer: thickness_nm must be positive');
end
if exist(material_file, 'file') ~= 2
    error('add_layer: material file not found: %s', material_file);
end

layer.material_file = material_file;
layer.thickness_nm  = thickness_nm;

% Derive a short label from the filename for CSV naming.
% Expects pattern like  n_Au_cxro.txt  ->  'Au'
[~, fname, ~] = fileparts(material_file);
parts = strsplit(fname, '_');
if numel(parts) >= 2
    layer.label = parts{2};
else
    layer.label = fname;
end

stack.layers{end+1}   = layer;
stack.total_height_nm = stack.total_height_nm + thickness_nm;

end
