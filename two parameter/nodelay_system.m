function states_dot = nodelay_system(t,states, physical_parameters, a_normalized,omega)  

m = physical_parameters(1); 
l = physical_parameters(2); 
kappa = physical_parameters(3); 
k1 = physical_parameters(4); 
k2 = physical_parameters(5);  
M = physical_parameters(6); 
C = physical_parameters(7);  
K = physical_parameters(8);  
g = physical_parameters(9);


a = a_normalized * sqrt(omega^2*C^2 + K^2);
timescale = omega / (4*pi); %to scale outputs to have a period of 1

%dydt = [thetadot; thetaddot; ydot; yddot]
%y = [theta; thetadot; y; ydot]

states_dot = zeros(4,1);
states_dot(1) = states(2);
states_dot(3) = states(4);

LHS_matrix = [m*l^2, m*l*sin(states(1)); m*l*sin(states(1)), M+m].*timescale^2;
RHS_matrix = [-kappa*states(2)*timescale - m*g*l*sin(states(1));
              -C*states(4)*timescale - K*states(3) - m*l*(states(2)^2)*cos(states(1))*(timescale^2) + a*cos(4*pi*t)];
          
states_dot([2,4]) = LHS_matrix \ RHS_matrix;
end