function [time,states] = extra(a_normalized,omega,Tf,points,physical_parameters) 
history = [pi/24; 0; 0; 0];
tspan = [0 Tf];
tau = physical_parameters(10);
sol = dde23(@ddex1de,tau,history,tspan,[],omega,a_normalized,physical_parameters);
time = linspace(0,Tf,points);
states = deval(sol,time);
end

function states_dot = ddex1de(t,states,Z,omega,a_normalized,physical_parameters)  

m = physical_parameters(1); 
l = physical_parameters(2); 
kappa = physical_parameters(3); 
k1 = physical_parameters(4); 
k2 = physical_parameters(5);  
M = physical_parameters(6); 
C = physical_parameters(7);  
K = physical_parameters(8);  
g = physical_parameters(9);
tau = physical_parameters(10);
a = a_normalized*sqrt(omega^2*C^2+K^2);

timescale = (omega/(4*pi)); %Timescale to scale inputs to period 1

states_dot = zeros(4,1);
states_dot(1) = states(2);
states_dot(3) = states(4);
ya = Z(4);

LHS_matrix = [m*l^2, m*l*sin(states(1)); m*l*sin(states(1)), M+m].*timescale^2;
RHS_matrix = [-kappa*states(2)*timescale - m*g*l*sin(states(1)) - m*l*ya*sin(states(1)*(timescale^2));
              -C*states(4)*timescale - K*states(3) - m*l*(states(2)^2)*cos(states(1))*(timescale^2) + a*cos(4*pi*t) - m*l*ya*sin(states(1)*(timescale^2))];

states_dot([2,4]) = LHS_matrix \ RHS_matrix;
end