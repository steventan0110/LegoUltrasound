
% clear all;
addpath k-Wave
rng('shuffle');

%% input data
% transducer
c0 = 1500;                  % medium speed of sound [m/s]
f0 = 1.0e6;                 % transducer central frequency [Hz]
% bw = 0.8;                   % pulse bandwidth
% p0 = 0.5e6;                 % pulse amplitude [Pa]
Nt = 96;                    % number of transmit elements
Nr = 96;                    % number of receive elements
lambda = c0/f0;             % wavelength [m]
% origin = [0 0 0];           % origin [m,m,m]

% domain
rho0 = 990;                 % medium density [kg/m3]
pml_x_size = 8;
pml_y_size = 4;
pml_z_size = 8;
Nx = 128 - 2 * pml_x_size;  % number of elements in the domain
Ny = 64 - 2 * pml_y_size;
Nz = 128 - 2 * pml_z_size;  % number of elements in the domain


%% spatial and temporal grid
% spatial grid
h = 4e-4;                       % desired spatial step [m]
x = (1:Nx)*h;                   % x-vector [m]
y = (1:Ny)*h;                   % y-vector [m]
z = (1:Nz)*h;                   % z-vector [m]
kgrid = kWaveGrid(Nz,h,Nx,h,Ny,h);    % k-wave structure

% temporal grid
t_end = sqrt((Nz + 2*pml_z_size)^2 + (Nx + 2*pml_x_size)^2)*h/c0 + 4/f0;             % maximum value of time vector [s]
CFL = 0.25;                     % CFL number
[kgrid.t_array, delta] = makeTime(kgrid, c0, CFL, t_end);
                                % k-wave structure
t = kgrid.t_array;              % t-vector (s)
Fs = 1/delta;                   % sampling frequency (Hz), Fs = c0/(CFL*h)

% stability check
assert(f0 < Fs/2, 'The excitation frequency exceeds Nyquist limit in time.');
assert(f0<c0/2/h,'The excitation frequency exceeds Nyquist limit in space.');


%% transducer definition
transducer.number_elements = Nr;                    % total number of transducer elements
transducer.element_width = round(0.2798e-3/h);      % width of each element [grid points]
transducer.element_length = round(4e-3/h);          % length of each element [grid points]
transducer.element_spacing = round(0.025e-3/h);     % spacing (kerf width) between the elements [grid points]
transducer.radius = inf;                            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
+ (transducer.number_elements - 1) * transducer.element_spacing;

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                        % sound speed [m/s]
% transducer.focus_distance = 20e-3;                % focus distance [m]
transducer.elevation_focus_distance = 16e-3;        % focus distance in the elevation plane [m]
transducer.steering_angle = 0; 
% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Nx/2-transducer_width/2, Ny/2 - transducer.element_length/2]);
% % append input signal used to drive the transducer
% transducer.input_signal = input_signal;
transducer.receive_apodization = 'Hanning';
% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);


%% source definition
% define properties of the input signal (only when transducer is the acoustic source)
source_strength = 1e6;          % [MPa]
tone_burst_freq = 0.5e6;        % [Hz]
tone_burst_cycles = 4;

source.p = source_strength .* toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);


%% data generation
num_iter = 1500;                      % number of samples to create
input_args = {'PMLSize',[pml_z_size,pml_x_size,pml_y_size],'PMLInside',...
    false,'PlotSim',false,'RecordMovie',false,'DataCast','gpuArray-single'};
    % channel_data = zeros(length(t), Nt, Nr);

for iter = 1485:num_iter
    % make directionary
    pathname = strcat(pwd, '/data/', num2str(iter));
    mkdir(pathname);
    % domain definition
    medium_mask = makeRandomLesions([Nz, Nx, Ny]);
    sos_num = max(medium_mask(:));
    sos_bg = normrnd(1600, 10);
    sos_prob = 1 - round(rand([1, sos_num])*0.6);
    sos_medium = sos_prob .* (1570 + 60 * (rand([1, sos_num]) - 0.5)) + ...
        (1 - sos_prob) .* (1625 + 30 * (rand([1, sos_num]) - 0.5));
    att_bg = normrnd(0.186, 0.014);
    att_medium = normrnd(0.142, 0.008, [1, sos_num]);

    % medium matrix
    medium.sound_speed = sos_bg * ones(Nz,Nx,Ny);   % sound speed [m/s]
    medium.density = rho0 * ones(Nz,Nx,Ny);         % density [kg/m3]
    medium.alpha_coeff = att_bg * ones(Nz,Nx,Ny);   % [dB/(MHz^y cm)]
    for i = 1:sos_num
        medium.sound_speed(medium_mask == i) = sos_medium(i);
        medium.density(medium_mask == i) = rho0 + 20*(rand-0.5);
        medium.alpha_coeff(medium_mask == i) = att_medium(i);
    end

    medium.alpha_power = 1.05;                  % non-linear relationship between attenuation and frequency
    medium.alpha_mode = 'no_dispersion';
    medium. BonA = 6;                           % non-linearity
    
    % save the labels for sos and attenuation
    medium_sos = medium.sound_speed;
    medium_att = medium.alpha_coeff;
    medium_lsos = medium.sound_speed(:, :, Ny/2);
    medium_latt = medium.alpha_coeff(:, :, Ny/2);
    save(strcat(pathname, '/sos.mat'), 'medium_sos');
    save(strcat(pathname, '/att.mat'), 'medium_att');
    save(strcat(pathname, '/lb_sos.mat'), 'medium_lsos');
    save(strcat(pathname, '/lb_att.mat'), 'medium_latt');

    % transmit signal element by element
    for f = 1:Nt
        source_mask = zeros(Nz, Nx, Ny);
        source_mask(Nz, Nx/2-transducer_width/2 + (f-1)*transducer.element_width, Ny/2) = 1;
        source.p_mask = source_mask;
        filename = strcat(pathname,'/',num2str(f),'.h5');
        sensor_data = kspaceFirstOrder3D(kgrid,medium,source,transducer,input_args{:},'SaveToDisk',filename);
    %    sensor_data = kspaceFirstOrder3DG(kgrid,medium,source,transducer,input_args{:});
    %    channel_data(:, :, f) = sensor_data';
    end
end
% save(date)
% figure, imagesc(channel_data(:, :, 1)), colormap(gray);


