% channel data with k-wave
clear;
addpath .\beamforming

%% input data

% domain
c0 = 1540;            % medium speed of sound [m/s]
rho0 = 1020;          % medium density [kg/m3]
Lx = 60e-3;           % width of domain (m)
Lz = 50e-3;           % depth of domain (m)

% transducer
f0 = 2.0e6;           % transducer central frequency [Hz]
bw = 0.8;             % pulse bandwidth
p0 = 0.5e6;           % pulse amplitude [Pa]
N = 128;              % number of elements
lambda = c0/f0;       % wavelengsth [m]
pitch = lambda/2;     % element pitch [m]
origin = [0 0 0];     % origin [m,m,m]
focus = [0 0  15e-3]; % transmit focus [m,m,m]


%% Spatial and temporal grid

% spatial grid
h = lambda/6;                  % desired spatial step (m)
Nx = round(Lx/h);              % number of elements in the domain
Nz = round(Lz/h);              % number of elements in the domain
x = ((0:(Nx-1))-(Nx-1)/2)*h;   % x-vector [m]
z = (0:(Nz-1))*h;              % z-vector [m]
[xx,zz] = meshgrid(x,z);       % x and z matrixes [m]
kgrid = makeGrid(Nz,h,Nx,h);   % k-Wave structure
PML_Size = round(4*lambda/h);  % size of the PML layers in grid points

% temporal grid
t_end = 2*(Lz/c0+3/f0);        % maximum value of time vector [s]
CFL = 0.3;                     % CFL number, CFL = c0*delta/h
[kgrid.t_array, delta] = makeTime(kgrid, c0, CFL, t_end); % k-Wave structure
t = kgrid.t_array;             % t-vector (s)
Fs = 1/delta;                  % sampling frequency (Hz)

% stability check
assert(f0<Fs/2,'The excitation frequency exceeds Nyquist limit in time.')
assert(f0<c0/2/h,'The excitation frequency exceeds Nyquist limit in space.')


%% domain definition
% wires
point_z = 25e-3;             % z-position of the scatterers [m]
point_x = 0;                 % x-position of the scatterers [m]

% snapping points to grid, this code is valid for more than one points
point_sample_z = (round((point_z-z(1))/h)+1); 
point_sample_z = point_sample_z(point_sample_z < Nz);
point_sample_x = (round((point_x-x(1))/h+1)); 
point_sample_x = point_sample_x(point_sample_x < Nx);
[wzwz, wxwx] = meshgrid(point_sample_z, point_sample_x);
point_mask = zeros(Nz,Nx);
for n = 1:length(wzwz(:))
    point_mask = point_mask + makeDisc(Nz,Nx,wzwz(n),wxwx(n),round(lambda/8/h));
end

% medium matrix
medium.sound_speed               = c0*ones(Nz, Nx);     % sound speed [m/s]
medium.density                   = rho0*ones(Nz, Nx);   % density [kg/m3]
medium.sound_speed(point_mask==1) = c0;                  % sound speed in wire [m/s]
medium.density(point_mask==1)     = 2*rho0;              % density of wire [kg/m3]


%% source definition
% probe geometry and pulse definition
geom=[pitch*((1:N)-(N+1)/2).' zeros(N,2)]; % probe's geometry
t0=-5/f0:delta:5/f0;                       % pulse time vector [s]
pulse=p0*gauspuls(t0,f0,bw);               % generated pulse [Pa]

% transmit focus & apodization
tx_apodization=hamming(N);                            % transmit apodization
tx_delay=(sqrt(sum((origin-focus).^2,2))-sqrt(sum((geom-ones(N,1)*focus).^2,2)))/c0;
t00=(min(tx_delay)-5/f0):delta:(max(tx_delay)+5/f0);  % pulse time vector [s]
t=t+min(t00);                                         % redefining time vector to compensate for pulse length;

% snapping transducer to grid and assigning temporal signals
pitch_n=round(pitch/h);            % quantified pitch
[~,n0]=min(abs(x+N*pitch/2));   % initial transducer position
source.p=[]; 
source.p_mask=zeros(Nz,Nx); 
for n=1:N
    nx=n0+1+pitch_n*(n-1)+(0:(pitch_n-1));
    [~,nz]=min(abs(z-geom(n,3)));
    source.p_mask(nz,nx)=1;
    source.p=[source.p; repmat(tx_apodization(n)*interp1(t0,pulse,t00-tx_delay(n),'linear',0),pitch_n,1)];
end


%% sensor definition
sensor.mask=source.p_mask;

%% Run the simulation
input_args = {'PlotSim', false, 'PMLSize', PML_Size, 'PMLInside', false, 'PlotFreq', 10};% 'RecordMovie', true};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
p_sensor_data = permute(sensor_data,[2 1]);

%% Convert sensor data to RF data
% converting sensir data to channel data
rf_channel_data = zeros(size(p_sensor_data,1),N);
for n = 1:N
    rf_channel_data(:,n) = rf_channel_data(:,n)+sum(p_sensor_data(:,pitch_n*(n-1)+(1:pitch_n)),2);
end

cd = rf_channel_data(abs(round((min(tx_delay)-5/f0)/delta)) + 2:end, :);
% showing rf-channel data
figure
imagesc(1:N,t,cd); colormap gray;
title('RF channel data')

%% beamforming
% beamformed_data = DAS_single_2(cd,N*3,pitch*1e3,c0,1);
[S, N] = size(cd);

post_data = zeros(S, N);
sampleSpacing = c0/Fs;
for j = 1:S
    for i = 1:N
        delay = round((sqrt(((i-64)*pitch)^2 + (j/2*sampleSpacing)^2) - j/2*sampleSpacing) / sampleSpacing);
        if delay + j < S
            post_data(j, i) = post_data(j, i) + cd(j + delay, i);
        end
    end
end
Aline = sum(post_data, 2);
figure, plot(Aline);
figure, imagesc(post_data);

% % post processing
% post_data = abs(hilbert(beamformed_data));
% post_data = post_data/max(max(post_data));
% % show image
% figure
% imagesc(db(post_data),[-20,0]);
% colormap(hot);

