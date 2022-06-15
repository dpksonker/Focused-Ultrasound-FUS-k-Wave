
clear
clc
clearvars;

% medium properties White Matter
c0 = 2652;       %% Maximum medium sound speed
% source parameters
source_f0       = 0.5e6;      % source frequency [Hz]
source_roc      = 47.16e-3;    % bowl radius of curvature [m]
source_diameter = 50e-3;    % bowl aperture diameter [m]
source_amp      = 0.07e6;      % source pressure [Pa]

% grid parameters
axial_size      = 80e-3;    % total grid size in the axial dimension [m]
lateral_size    = 60e-3;    % total grid size in the lateral dimension [m]

% computational parameters
ppw             = 11;        % number of points per wavelength
record_periods  = 1;        % number of periods to record
cfl             = 0.1;      % CFL number
source_x_offset = 4;       % grid points to offset the source

% =========================================================================
% GRID
% =========================================================================

% calculate the grid spacing based on the PPW and F0

dx =   c0/ (ppw * source_f0);   % [m]

% compute the size of the grid

Nx = roundEven(lateral_size / dx);
Ny = roundEven(lateral_size / dx);
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
medium.alpha_power = 1.1;    %  Power law exponent

for ii = 1:Nx

    for jj = 1:Ny

        for kk = 1:Nz

            %                 rr = ceil(sqrt((Nx/2 - ii)^2 + (Ny/2 - jj)^2) - 20);
            rr = sqrt((ii - (Nx/2))^2 + (jj - (Ny/2))^2 + (kk - (Nz/2+ 26))^2);

            if (60 < rr) && (rr <= 66)
                %                scalp
                medium.sound_speed(ii,jj,kk) = 1537;	% [m/s]
                medium.density(ii,jj,kk) = 1116;       % [kg/m^3]
                medium.alpha_coeff(ii,jj,kk) = 0.85;  % [dB/(MHz^y cm)] Power Law Prefactor

            elseif (46 < rr) && (rr <= 60)
                %                skull
                medium.sound_speed(ii,jj,kk) = 2652;	% [m/s]
                medium.density(ii,jj,kk) = 1796;       % [kg/m^3]
                medium.alpha_coeff(ii,jj,kk) = 4.11;  % [dB/(MHz^y cm)] Power Law Prefactor

                %                     20 dB
            elseif (40 < rr) && (rr <= 46)
                %                     csf
                medium.sound_speed(ii,jj,kk) = 1500;	% [m/s]
                medium.density(ii,jj,kk) = 1000;       % [kg/m^3]
                medium.alpha_coeff(ii,jj,kk) = 0.002;  % [dB/(MHz^y cm)]  Power Law Prefactor
                % 0.75 is the      % % absorption  coefficient

            elseif rr <= 40
                %                     average brain
                medium.sound_speed(ii,jj,kk) = 1500;	% [m/s]
                medium.density(ii,jj,kk) = 1000;       % [kg/m^3]
                medium.alpha_coeff(ii,jj,kk) = 0.10;  % [dB/(MHz^y cm)]  Power Law Prefactor
                % 0.75 is the      % % absorption  coefficient


            end


        end
    end
end
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


% Convert the medium acoustic properties from dB to Neper
alpha_np = db2neper(medium.alpha_coeff, medium.alpha_power) * ...             
    (2 * pi * source_f0).^medium.alpha_power;


% % reshape data for pressure

p = reshape(p, Nx, Ny, Nz);

% Crop the pressure data for Acousto-optics

% n_p = p(11:115,11:115,43:147);
% phase = reshape(phase, Nx, Ny, Nz);
% n_ph = phase(11:115,11:115,43:147);
% u_x = reshape(sensor_data.ux_rms,Nx,Ny,Nz);
% u_xx = u_x(11:115,11:115,43:147);
% u_y = reshape(sensor_data.uy_rms,Nx,Ny,Nz);
% u_yy = u_x(11:115,11:115,43:147);
% u_z = reshape(sensor_data.uz_rms,Nx,Ny,Nz);
% u_zz = u_x(11:115,11:115,43:147);

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


% % plot the acoustic pressure
figure(20)
% subplot(3, 3, 1);
imagesc(kgrid.z_vec * 1e3, kgrid.x_vec * 1e3, squeeze(p(:,62,:)*1e-6));
h = colorbar;
xlabel(h, '[MPa]');
ylabel('x (mm)');
xlabel('z (mm)');
axis image;

% % set colormap and enlarge figure window
colormap(hot(256));
scaleFig(1.5, 1);
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',45,'FontWeight','Bold')

h = colorbar('Ticks',[ 0.1 0.2 0.3 0.4 0.5],'TickLabels',{'0.1','0.2','0.3','0.4','0.5'})
ylabel(h,'P_0 (MPa)','FontSize',45);
% %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
