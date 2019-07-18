%% channel data with k-wave
% 2D k-wave simulation with an array of 128 source immitating transducer's
% imgae. Each one of them emits wave generated toneBurst function which
% mimics the behavior of real ultrasound probe used in the medical
% imaging.
close all
clear;
clc;
addpath k-Wave %the k-Wave is stored directly in LegoUltrasound folder, otherwise needs change
addpath medium_def %for the medium definition
%% input data

% domain
c0 = 1540;            % medium speed of sound (in body tissue) [m/s]
rho0 = 1020;          % medium density [kg/m3]
Lx = 70e-3;           % width of domain (m)
Lz = 60e-3;           % depth of domain (m)

% transducer
f0 = 2.0e6;           % transducer central frequency [Hz]
bw = 0.8;             % pulse bandwidth
p0 = 0.5e6;           % pulse amplitude [Pa]
N = 128/2;            % number of elements， use 128 cause so many simulation to run
lambda = c0/f0;       % wavelengsth [m]
pitch = lambda/2*2;     % element pitch [m]
origin = [0 0 0];     % origin [m,m,m]
focus = [0 0  15e-3]; % transmit focus [m,m,m]


%% Spatial and temporal grid

% spatial grid
h = lambda/6;                  % desired spatial step (m)
Nx = round(Lx/h);              % number of elements in the domain
Nz = round(Lz/h);              % number of elements in the domain
x = ((0:(Nx-1))-(Nx-1)/2)*h;   % x-vector [m] since the center is in the middle
z = (0:(Nz-1))*h;              % z-vector [m]
[xx,zz] = meshgrid(x,z);       % x and z matrixes [m]
kgrid = makeGrid(Nz,h,Nx,h);   % k-Wave structure

PML_Size = round(4*lambda/h)+1;  % size of the PML layers in grid points

% temporal grid
diagonal = (Lz^2 + Lx^2)^0.5;   %the length of the diagonal in the grid
t_end = 2*(diagonal/c0);  % maximum value of time vector [s]
CFL = 0.6;                     % CFL number, CFL = c0*delta/h
[kgrid.t_array, delta] = makeTime(kgrid, c0, CFL, t_end); % k-Wave structure
t = kgrid.t_array;             % t-vector (s)
Fs = 1/delta;                  % sampling frequency (Hz)


% stability check
assert(f0<Fs/2,'The excitation frequency exceeds Nyquist limit in time.')
assert(f0<c0/2/h,'The excitation frequency exceeds Nyquist limit in space.')


%% domain definition
% This part defines the medium, which will be replaced by medium retrived
% from get_medium function. For now use a point
[image, ind] = generate(Nz, Nx);
point_mask = zeros(Nz,Nx);
point_mask(image ~= 0) = 1;
%imshow(point_mask);
% medium matrix
medium.sound_speed               = c0*ones(Nz, Nx);     % sound speed [m/s]
medium.density                   = rho0*ones(Nz, Nx);   % density [kg/m3]
medium.sound_speed(point_mask==1) = c0;                  % sound speed in wire [m/s]
medium.density(point_mask==1)     = 2*rho0;              % density of wire [kg/m3]


%% source definition and Beamforming for Bmode image process:

pitch_n=round(pitch/h);         % quantified pitch
[~,n0]=min(abs(x+N*pitch/2));   % initial transducer position

% probe geometry and pulse definition
geom=[pitch*((1:N)-(N+1)/2).' zeros(N,2)]; % probe's geometry
t0=-5/f0:delta:5/f0;                        % pulse time vector [s]
pulse=p0*gauspuls(t0,f0,bw);               % generated pulse [Pa]

numPics = 10; % define the number of images I want
for pic = 1: numPics
    Bscan = zeros(1400, N); %Bscan compose of N number of Aline
    
    for i = 1:N  %iterate through the elements for transmission and beamform
        %% source definition
        %redefine the origin and focus:
        nx= x(n0) + (i-1)*pitch;
        origin = [nx 0 0];     % origin [m,m,m]
        focus = [nx 0 15e-3]; % transmit focus [m,m,m], focus depth is never changed
        
        % transmit focus & apodization
        tx_apodization=hamming(N);                             % transmit apodization
        d_vert =  sqrt(sum((origin-focus).^2,2)); %the vertical distance, which is constant
        
        d_slope = sqrt(sum((geom-ones(N,1)*focus).^2,2)); %the longest side of the triangle for each transmission element
        tx_delay=(d_vert-d_slope)/c0;
        t00=(min(tx_delay)-5/f0):delta:(max(tx_delay)+5/f0);  % pulse time vector [s]
        t = kgrid.t_array;
        t=t+min(t00);                                         % redefining time vector to compensate for pulse length;
        
        % snapping transducer to grid and assigning temporal signals
        source.p=[];
        source.p_mask=zeros(Nz,Nx);
        
        for n=1:N %define the source mask for simulation
            nx_s=n0+1+pitch_n*(n-1)+(0:(pitch_n-1));
            [~,nz]=min(abs(z-geom(n,3)));
            source.p_mask(nz,nx_s)=1;
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
        rf_channel_data = zeros(size(p_sensor_data,1),N); %N channels since N is # of source
        for n = 1:N
            rf_channel_data(:,n) = rf_channel_data(:,n)+sum(p_sensor_data(:,pitch_n*(n-1)+(1:pitch_n)),2);
        end
        
        cd = rf_channel_data(abs(round((min(tx_delay)-5/f0)/delta)) + 2:end, :);
        % showing rf-channel data
        %     figure
        %     imagesc(1:N,t,cd); colormap gray;
        %     title('RF channel data')
        
        %% beamforming
        [S, N_b] = size(cd);
        post_data = zeros(S, N_b);
        sampleSpacing = c0/Fs;
        
        for j = 1:S
            for k = 1:N_b
                halfapt = max(max(i-1, N-i));
                % apodization
                win = hamming(halfapt*2 + 1);
                %distance from reference point to Nth element:
                d1 = (i-k)*pitch; %the horizontal distance
                d2 = j/2*sampleSpacing; %the vertical distance, divide by 2 cause it's a round trip
                dif = (sqrt(d1^2+d2^2) - d2);
                delay = round(dif/sampleSpacing);
                if delay + j < S
                    post_data(j, k) = post_data(j, k) + win(k - i + halfapt + 1) * cd(j + delay, k);
                end
            end
        end
        Aline = sum(post_data, 2); %Aline's dimension would be [S, 1]
        if size(Aline,1) < 1400
            Aline(1400,1) = 0;
        end
        Bscan(:, i) = Aline(1:1400,:); %concatenate Aline into Bscan image
        
        
    end
    
%     output the image for testing
%     figure
%     imagesc(logCompression(abs(hilbert(Bscan)), 0.01)), colormap(gray);
pathname = strcat(pwd, '/data/', num2str(pic));
mkdir(pathname);
save(strcat(pathname, 'Bscan.mat'),'Bscan');
end

% figure, plot(Aline);
% figure, imagesc(post_data);

% % post processing
% post_data = abs(hilbert(beamformed_data));
% post_data = post_data/max(max(post_data));
% % show image
% figure
% imagesc(db(post_data),[-20,0]);
% colormap(hot);


