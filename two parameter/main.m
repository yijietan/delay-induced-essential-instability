% Parameters for numerical simulation
runtime = 20;
points = 2000;
modes = 2;
stepsize = 1e-4;
tolerance = 1e-4;
continuation_stepsize = 0.01;
continuation_points = 50;              
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
r = 0.1;
physical_parameters = [m; l; kappa; k1; k2; M; C; K; g; tau; r];

% Initial guess parameters
omega = 17; 
a_normalized = 0.01; 
initial_conditions = [pi/24; 0; 0; 0];          % [theta; theta_dot; y; y_dot]
time_array = linspace(0,runtime,points);

% returns initial guess of the solution based on outputs of system w/o delay
z_initialguess = initial_guess(runtime, physical_parameters, omega, a_normalized, modes, initial_conditions, time_array);

% returns the first point on the continuation curve with a Newton iteration
z_firstpoint = Newton_firstpoint(z_initialguess, a_normalized, omega, numerical_parameters, physical_parameters, initial_conditions, time_array);

% extends z vector with omega to perform pseudo-arclength continuation
z_0 = [z_firstpoint; omega];

% performs Newton continuation with pseudo-arclength continuation
z_new = Newton_continuation(z_0, numerical_parameters, physical_parameters, initial_conditions, time_array);

%% Plot continuation curve
f1 = figure;
plot(z_new(12,:),z_new(11,:));
grid on;
xlim([12 17]);
title('Continuation curve of period doubling');
xlabel('\Omega (rad/s)');
ylabel('a/(\Omega^2C^2 + K^2)^{1/2} (m)');
saveas(f1,'continuation_curve_test','png');