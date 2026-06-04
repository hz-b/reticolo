clear; warning('off', 'all');

base    = fileparts(mfilename('fullpath'));
oc_path = fullfile(base, '..', '..', 'Optical_Constants');
addpath(fullfile(base, '..', '..', 'helpers'));


% Geometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_res_nm = 0.5;
z_res_nm = 0.5;

grPeriod_lpermm      = 600;
grBlazeAngle_deg     = 0.73;
grAntiBlazeAngle_deg = 5.60;

period_nm  = 1e6 / grPeriod_lpermm;
tan_blaze  = tand(grBlazeAngle_deg);
tan_anti   = tand(grAntiBlazeAngle_deg);
w_blaze_nm = period_nm / (1 + tan_blaze / tan_anti);
depth_nm   = w_blaze_nm * tan_blaze;

fprintf('Derived groove depth: %.4f nm\n', depth_nm);

grating = build_grating('blazed', period_nm, depth_nm, ...
                         grBlazeAngle_deg, grAntiBlazeAngle_deg, x_res_nm);


% Layer stack %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

substrate_file = fullfile(oc_path, 'n_Si_cxro.txt');

stack = build_stack(grating);
stack = add_layer(stack, fullfile(oc_path, 'n_Au_cxro.txt'), 31);   % Au coating
% stack = add_layer(stack, fullfile(oc_path, 'n_C_cxro.txt'), .8);   % C contamination


% Sweep parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sweep.type      = 'energy';
sweep.values    = 50:5:1000;
% sweep.alpha_deg = 4;        % fixed grazing incidence angle in degrees
sweep.Cff     = 1.5;     % to use Cff-based angle instead


% Solver options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.FourierOrders = 11;
options.pol           = -1;          % -1 = TM,  +1 = TE
options.GR_Order      = -1;
options.z_res_nm      = z_res_nm;
options.reticolo_path = fullfile(base, '..', '..', '..', 'V9', 'reticolo_allege_v9');
options.output_dir    = fullfile(base, 'Results');
options.oc_path       = oc_path;
options.verbose       = true;


% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = run_rcwa(stack, substrate_file, sweep, options);


% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_results(stack, substrate_file, sweep, options, results);
