% RUN_RCWA  Run a RETICOLO RCWA sweep over photon energy or incidence angle.
%
% USAGE
%   results = run_rcwa(stack, substrate_file, sweep, options)
%
% INPUTS
%   stack           struct from build_stack / add_layer
%
%   substrate_file  path to CXRO optical constants file for the substrate
%
%   sweep           struct describing what to iterate over:
%
%     Energy sweep with Cff mode (alpha computed per energy):
%       sweep.type    = 'energy'
%       sweep.values  = [50 : 5 : 1000]   % eV
%       sweep.Cff     = 2.25
%
%     Energy sweep with fixed alpha:
%       sweep.type      = 'energy'
%       sweep.values    = [50 : 5 : 1000]
%       sweep.alpha_deg = 1.5              % grazing incidence angle
%
%     Angle sweep at fixed energy:
%       sweep.type        = 'alpha'
%       sweep.values      = [0.5 : 0.1 : 5.0]   % grazing incidence angles (deg)
%       sweep.energy_eV   = 500                  % fixed photon energy
%
%   options         struct (all fields optional, defaults shown):
%       .FourierOrders   11
%       .pol             -1     % -1 = TM, +1 = TE
%       .GR_Order        -1     % diffraction order to extract
%       .z_res_nm        0.5
%       .reticolo_path   fullfile(pwd, 'V9', 'reticolo_allege_v9')
%       .output_dir      pwd
%       .verbose         true
%
% OUTPUT
%   results   struct with fields:
%       .sweep_values       the iterated variable (eV or deg)
%       .sweep_label        'PhotonEnergy_eV' | 'GrazingAngle_deg'
%       .efficiency         diffraction efficiency for GR_Order
%       .alpha_deg          grazing incidence angle used at each point
%       .beta_deg           exit grazing angle for GR_Order
%       .csv_file           path of the written CSV