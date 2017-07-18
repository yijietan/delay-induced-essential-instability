%Parameters for numerical simulation
runtime = 20;
points = 2000;
modes = 6;
stepsize = 1e-4;
tolerance = 1e-4/2;
continuation_stepsize = 1e-3;
continuation_points = 349;              
direction = -1;                                 % -1 for leftwards, 1 for rightwards
numerical_parameters = [runtime; points; modes; stepsize; tolerance; continuation_stepsize; continuation_points; direction]; 

% Physical parameters
m = 0.27; 
l = 0.1955; 
kappa = 7.5e-3;
k1 = 4; 
k2 = 4; 
M = 0.2; 
C = 0.1; 
K = 200; 
g = 9.81;
tau = 0.01;
r = 0;
physical_parameters = [m; l; kappa; k1; k2; M; C; K; g; tau; r];

% Initial guess parameters
omega = 13; 
a_normalized = 0.014; 
initial_conditions = [pi/24; 0; 0; 0];          % [theta; theta_dot; y; y_dot]
time_array = linspace(0,runtime,points);

% returns initial guess of the solution based on outputs of system w/o delay
z_initialguess = initial_guess(runtime, physical_parameters, omega, a_normalized, modes, initial_conditions, time_array);

% returns the first point on the continuation curve with a Newton iteration
z_firstpoint = Newton_firstpoint(z_initialguess, a_normalized, omega, numerical_parameters, physical_parameters);

[maxtheta, z_new] = Newton_continuation(z_firstpoint, numerical_parameters, physical_parameters,omega,time_array);

%% plot the continuation curve
f1 = figure;
plot(z_new(end,:), maxtheta);
grid on;
xlabel('a/(\Omega^2C^2 + K^2)^{1/2} (m)');
ylabel('max \theta (rad)');
title('One parameter continuation for \Omega = 13 rad/s');
saveas(f1,'one parameter continuation', 'png');