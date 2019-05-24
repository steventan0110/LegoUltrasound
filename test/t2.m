% define the physical properties of the transducer
tr.number_elements = 128; % total number of transducer elements
tr.element_width = 1; % width of each element [grid points]
tr.element_length = 12; % length of each element [grid points]
tr.element_spacing = 0; % spacing between the elements [grid points]

tr.radius = Inf; % radius of curvature of the transducer [m]
% define the position of the transducer [grid points]
tr.position = [1, 20, 20];
% define the properties used to derive the beamforming delays
tr.sound_speed = 1540; % sound speed [m/s]
tr.focus_distance = 20e-3; % focus distance [m]
tr.steering_angle = 10; % steering angle [degrees]
% define the apodization
tr.transmit_apodization = 'Rectangular';
tr.receive_apodization = 'Hanning';
% define the transducer elements that are currently active
tr.active_elements = zeros(tr.number_elements, 1);
tr.active_elements(1:32) = 1;
% define the input signal used to drive the transducer
tr.input_signal = input_signal;
% create the transducer using the defined settings
transducer = makeTransducer(kgrid, tr);
% display the transducer using a 3D voxel plot
transducer.plot;
% print a list of transducer properties to the command line
transducer.properties;

% run a simulation using the same transducer as both source and sensor
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer);
% form the recorded sensor data in a single scan line based on the current
% beamforming and apodization settings
scan_line = transducer.scan_line(sensor_data);