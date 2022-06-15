
clear  
clc
clearvars;

% medium properties White Matter
c0 = 1500;       %% Maximum medium sound speed
% source parameters
source_f0       = 0.5e6;      % source frequency [Hz]
source_roc      = 47.16e-3;    % bowl radius of curvature [m]
source_diameter = 50e-3;    % bowl aperture diameter [m]
source_amp      = 0.035e6;      % source pressure [Pa]

% grid parameters
axial_size      = 80e-3;    % total grid size in the axial dimension [m]
lateral_size    = 60e-3;    % total grid size in the lateral dimension [m]

% computational parameters
ppw             = 6;        % number of points per wavelength
record_periods  = 1;        % number of periods to record
cfl             = 0.1;      % CFL number
source_x_offset = 1;       % grid points to offset the source

% =========================================================================
% GRID
% =========================================================================

% calculate the grid spacing based on the PPW and F0

dx =  c0/ (ppw * source_f0);   % [m]

% compute the size of the grid

Ny = roundEven(lateral_size / dx);
Nx = roundEven(lateral_size / dx);
Nz = roundEven(axial_size / dx);

% create the computational grid

kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% compute points per temporal period

PPP = round(ppw / cfl);

% compute corresponding time spacing

dt = 1 / (PPP * source_f0);

% create the time array using an integer number of points per period

t_end = sqrt( kgrid.x_size.^2 + kgrid.y_size.^2 + kgrid.z_size.^2 )...
    / c0;     % total compute time [s] (this must be long enough to reach steady state)
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(c0 / (dx * source_f0))]);
disp(['CFL = ' num2str(c0 * dt / dx)]);

% =========================================================================
% SOURCE
% =========================================================================

% create time varying source
source_sig = createCWSignals(kgrid.t_array, source_f0, source_amp, 0);

% set bowl position and orientation
bowl_pos = [0, 0, kgrid.z_vec(1) + source_x_offset * kgrid.dx];
focus_pos = [0, 0, kgrid.z_vec(end)];

% create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.01, 'UpsamplingRate', 10);

% add bowl shaped element
karray.addBowlElement(bowl_pos, source_roc, source_diameter, focus_pos);
    
% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
    
% =========================================================================
% MEDIUM
% =========================================================================
% 
% define the properties of the propagation medium
                    
medium.sound_speed = 1500* ones(Nx, Ny,Nz); 	% [m/s]
medium.density = 1000* ones(Nx, Ny,Nz);        % [kg/m^3]
medium.alpha_coeff = 0.002* ones(Nx, Ny,Nz);  % [dB/(MHz^y cm)]  Power Law Prefactor
% medium.alpha_coeff = 0.21* ones(Nx, Ny,Nz);  % [dB/(MHz^y cm)]  Power Law Prefactor
medium.alpha_power = 1.1;    %  Power law exponent
% 

% --------------------
% SENSOR
% --------------------

% set sensor mask to record central plane, not including the source point
sensor.mask = ones(Nx, Ny, Nz);
% sensor.mask(source_x_offset + 2:end, :, Nz/2 + 1) = 1;

% record the pressure
sensor.record = {'p','u_rms'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * PPP + 1;

% =========================================================================
% SIMULATION
% =========================================================================

% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off', 'PlotScale', 'auto'};

sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:},'DataCast','single');

% =========================================================================
% % CALCULATE HEATING
% =========================================================================
% 
% extract amplitude from the sensor data
[p, phase] = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0);

alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...
    (2 * pi * source_f0).^medium.alpha_power;

% % % reshape data for pressure

p = reshape(p, Nx, Ny, Nz);

% Crop the pressure data for Acousto-optics

% n_p = p(11:111,11:111,60:160);
% phase = reshape(phase, Nx, Ny, Nz);
% n_ph = phase(11:111,11:111,60:160);
% u_x = reshape(sensor_data.ux_rms,Nx,Ny,Nz);
% u_xx = u_x(11:111,11:111,60:160);
% u_y = reshape(sensor_data.uy_rms,Nx,Ny,Nz);
% u_yy = u_x(11:111,11:111,60:160);
% u_z = reshape(sensor_data.uz_rms,Nx,Ny,Nz);
% u_zz = u_x(11:111,11:111,60:160);

Intensity = 0.5*p.^2 ./ (medium.density .* medium.sound_speed);

% filename = 'u_x.mat';
% save(filename,"u_xx",'-mat')
% filename = 'u_y.mat';
% save(filename,"u_yy",'-mat')
% filename = 'u_z.mat';
% save(filename,"u_zz",'-mat')
% filename = 'Pressure.mat';
% save(filename,"n_p",'-mat')
% filename = 'Phase.mat';
% save(filename,"n_ph",'-mat')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
figure(99)
% % plot the acoustic pressure
% subplot(3, 3, 1);
imagesc(kgrid.z_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p(:,(kgrid.Ny)/2+1,:)*1e-6));
h = colorbar;
xlabel(h, '[MPa]');
ylabel('x-position [mm]');
xlabel('z-position [mm]');
axis image;
title('Acoustic Pressure Amplitude water in Axial 0.5 MHz and 0.05 MPa');

% % set colormap and enlarge figure window
colormap(hot(256));
scaleFig(1.5, 1);

figure(990)
imagesc(squeeze(n_p(:,50,:)*1e-6));





