%2D grid with a sourse and sensor as the same line:
%the medium parameter is defined with ordinary value in water
%and bone - 2100

%need to define kgrid, medium, source, and sensor:
% create the computational grid
%define the grid as 256 * 256 
Nx = 256;
Ny = 256;
dx = 1;
dy = 1;
kgrid = makeGrid(Nx, dx, Ny, dy);
% number of grid points in the x (row) direction  
% number of grid points in the y (column) direction
% grid point spacing in the x direction [m]
% grid point spacing in the y direction [m]
% define the medium properties
medium.sound_speed = 1500*ones(Nx, Ny); % [m/s]
%half of the material is 2100m/s
medium.sound_speed(128:256,:) = 2100;     % [m/s]
medium.density = 1040;                  % [kg/m^3]
% define an initial pressure using makeLine
line_start_x = round(Nx/3);
line_end_x = round(2*Nx/3);
line_y = Ny/2;
line_mag = 3;
source.p0 = line_mag*makeLine(Nx, Ny, [line_start_x, line_y], [line_end_x, line_y]);

%define the binary mask for sensor:
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, :) = 1;

imagesc(double(sensor.mask));
% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PlotSim', true);

figure;
plot(sensor_data(1, :), 'b-');
plot(sensor_data(10,:), 'g-');
axis tight;
ylabel('Pressure');
xlabel('Time Step');
legend('Sensor Position 1');