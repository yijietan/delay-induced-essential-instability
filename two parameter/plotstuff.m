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
points = 1000;
Tf = 20;
omega = 20; 
a_normalized = 0.1; 
time_array = linspace(0,Tf,points);
initial_conditions = [pi/24; 0; 0; 0];   


% [~,states] = ode45(@(t, states)nodelay_system(t, states, physical_parameters, a_normalized, omega), time_array, initial_conditions);
% states = states.';
% plot(time_array,states(3,:))
% title('Plot of displacement of mass');
% xlabel('time (s)');
% ylabel('y(t) (m)');

%%
[time,states] = extra(a_normalized,omega,Tf,points,physical_parameters);

plot(time,states(3,:))
title('Plot of displacement of mass showing instability');
xlabel('time (s)');
ylabel('y(t) (m)');