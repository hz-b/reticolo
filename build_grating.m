function grating = build_grating(type, period_nm, depth_nm, varargin)
% BUILD_GRATING  Create a grating geometry descriptor.
%
% USAGE
%   grating = build_grating('blazed',      period_nm, depth_nm, ...
%                            blaze_deg, antiblaze_deg, x_res_nm)
%
%   grating = build_grating('trapezoidal', period_nm, depth_nm, ...
%                            duty_cycle, x_res_nm)
%
% OUTPUT  grating struct fields:
%   .type          'blazed' | 'trapezoidal'
%   .period_nm
%   .depth_nm
%   .x             [1 x Nx]  x coordinates across one period (nm)
%   .z_surface     [1 x Nx]  groove surface height at each x (nm)
%                            z=0 is the groove valley (substrate surface),
%                            z=depth_nm is the tip/top of the groove.
%   .params        struct carrying the original inputs for CSV naming

grating.type      = type;
grating.period_nm = period_nm;
grating.depth_nm  = depth_nm;

switch lower(type)

    % ------------------------------------------------------------------
    case 'blazed'
    % ------------------------------------------------------------------
    % varargin: blaze_deg, antiblaze_deg, x_res_nm
    if numel(varargin) < 3
        error('build_grating:blazed needs blaze_deg, antiblaze_deg, x_res_nm');
    end
    blaze_deg     = varargin{1};
    antiblaze_deg = varargin{2};
    x_res_nm      = varargin{3};

    tan_b = tand(blaze_deg);
    tan_a = tand(antiblaze_deg);

    % Blaze-face width from consistent two-angle geometry:
    %   depth = w_b * tan_b  and  depth = w_a * tan_a
    %   w_b + w_a = period
    w_blaze_nm = period_nm / (1 + tan_b / tan_a);
    x_peak_nm  = w_blaze_nm;

    % Guard against numerical edge cases
    edge = max(x_res_nm, period_nm * 1e-6);
    x_peak_nm = min(x_peak_nm, period_nm - edge);

    Nx = round(period_nm / x_res_nm) + 1;
    x  = linspace(0, period_nm, Nx);

    % Linear rise (blaze face) then hard drop (anti-blaze face)
    ctrl_x = [0,          x_peak_nm, period_nm];
    ctrl_z = [0,          depth_nm,  0         ];
    z_surface = interp1(ctrl_x, ctrl_z, x, 'linear');

    grating.params.blaze_deg     = blaze_deg;
    grating.params.antiblaze_deg = antiblaze_deg;
    grating.params.x_res_nm      = x_res_nm;

    % ------------------------------------------------------------------
    case 'trapezoidal'

    angle_deg   = varargin{1};
    width_ratio = varargin{2};
    x_res_nm    = varargin{3};

    W = width_ratio * period_nm;
    foot = depth_nm / tand(angle_deg);

    % --- validity check ---
    if 2*foot > (period_nm - W)
        error('Invalid trapezoid: sidewalls overlap');
    end

    % --- key points ---
    A = [W/2, 0];
    B = [W/2 + foot, depth_nm];
    C = [period_nm - (W/2 + foot), depth_nm];
    D = [period_nm - W/2, 0];

    Prf = [0,0; A; B; C; D; period_nm,0];

    % --- x grid (periodic-safe) ---
    Nx = round(period_nm / x_res_nm);
    x  = linspace(0, period_nm, Nx+1);
    x(end) = [];   % remove duplicate endpoint

    % --- interpolate safely ---
    z_surface = interp1(Prf(:,1), Prf(:,2), x, 'linear', 0);

    % --- store ---
    grating.x = x;
    grating.z_surface = z_surface;

    grating.params.angle_deg   = angle_deg;
    grating.params.width_ratio = width_ratio;
    grating.params.x_res_nm    = x_res_nm;

    otherwise
        error('Unknown grating type: %s.  Use ''blazed'' or ''trapezoidal''.', type);
end

grating.x         = x;
grating.z_surface = z_surface;

end