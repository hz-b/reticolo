function grating = build_grating(type, period_nm, depth_nm, varargin)
% BUILD_GRATING  Create a grating geometry descriptor.
%
% USAGE
%   grating = build_grating('blazed',      period_nm, depth_nm, ...
%                            blaze_deg, antiblaze_deg, x_res_nm)
%
%   grating = build_grating('trapezoidal', period_nm, depth_nm, ...
%                            width_ratio, angle_deg, x_res_nm)
%
%   width_ratio   groove BOTTOM width / period   (e.g. 0.67 = grWidthToD)
%   angle_deg     sidewall angle FROM HORIZONTAL (e.g. 75  = grTrapezoidAng_deg)
%   blaze_deg     blaze angle FROM HORIZONTAL (e.g. 0.729 = grBlazeAngle_deg)

grating.type      = type
grating.period_nm = period_nm
grating.depth_nm  = depth_nm

switch lower(type)

    % ------------------------------------------------------------------
    case 'blazed'
    % ------------------------------------------------------------------
    if numel(varargin) < 3
        error('build_grating: blazed needs blaze_deg, antiblaze_deg, x_res_nm');
    end
    blaze_deg     = varargin{1};
    antiblaze_deg = varargin{2};
    x_res_nm      = varargin{3};

    tan_b = tand(blaze_deg);
    tan_a = tand(antiblaze_deg);

    w_blaze_nm = period_nm / (1 + tan_b / tan_a);
    x_peak_nm  = w_blaze_nm;

    edge      = max(x_res_nm, period_nm * 1e-6);
    x_peak_nm = min(x_peak_nm, period_nm - edge);

    Nx = round(period_nm / x_res_nm) + 1;
    x  = linspace(0, period_nm, Nx);

    ctrl_x    = [0,        x_peak_nm, period_nm];
    ctrl_z    = [0,        depth_nm,  0         ];
    z_surface = interp1(ctrl_x, ctrl_z, x, 'linear');

    grating.params.blaze_deg     = blaze_deg;
    grating.params.antiblaze_deg = antiblaze_deg;
    grating.params.x_res_nm      = x_res_nm;

    % ------------------------------------------------------------------
    case 'trapezoidal'
    % ------------------------------------------------------------------
    % argument order: width_ratio, angle_deg, x_res_nm
    if numel(varargin) < 3
        error('build_grating: trapezoidal needs width_ratio, angle_deg, x_res_nm');
    end
    width_ratio = varargin{1}   % groove bottom width / period  (e.g. 0.67)
    angle_deg   = varargin{2}   % sidewall angle from horizontal (e.g. 75)
    x_res_nm    = varargin{3}

    if width_ratio <= 0 || width_ratio >= 1
        error('build_grating: width_ratio must be in (0,1)');
    end
    if angle_deg <= 0 || angle_deg > 90
        error('build_grating: angle_deg must be in (0,90]');
    end

    groove_floor_nm = width_ratio * period_nm;
    wall_foot_nm    = depth_nm / tand(angle_deg);   % horizontal component of each wall

    if 2 * wall_foot_nm > (period_nm - groove_floor_nm)
        error(['build_grating: trapezoid walls overlap. ' ...
               'Reduce depth, increase angle_deg, or reduce width_ratio.']);
    end

    % Land width at the top surface
    land_nm = period_nm - groove_floor_nm - 2 * wall_foot_nm;


    x_land_r  = land_nm / 2;
    x_floor_l = x_land_r  + wall_foot_nm;
    x_floor_r = x_floor_l + groove_floor_nm;
    x_land_l  = x_floor_r + wall_foot_nm;

    % Control points for interp1 (piecewise linear)
    ctrl_x = [0,         x_land_r,  x_floor_l, x_floor_r, x_land_l,  period_nm];
    ctrl_z = [depth_nm,  depth_nm,  0,         0,         depth_nm,  depth_nm ];

    Nx = round(period_nm / x_res_nm) + 1;
    x  = linspace(0, period_nm, Nx);

    z_surface = interp1(ctrl_x, ctrl_z, x, 'linear');

    grating.params.width_ratio = width_ratio;
    grating.params.angle_deg   = angle_deg;
    grating.params.x_res_nm    = x_res_nm;

    % ------------------------------------------------------------------
    otherwise
        error('build_grating: unknown type. Use ''blazed'' or ''trapezoidal''.');
end

grating.x         = x;
grating.z_surface = z_surface;

end