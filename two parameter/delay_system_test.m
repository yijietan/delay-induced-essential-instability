function states_dot = delay_system_test(t,states, physical_parameters, a_normalized,omega, z_in, modes) 

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

a = a_normalized * sqrt(omega^2*C^2 + K^2);
timescale = omega / (4*pi); %to scale outputs to have a period of 1

Fourier_sum = zeros(2*modes+2); 
Fourier_sum(1) = z_in(1);  
Fourier_sum(modes+2) = z_in(2*modes+2);

for i = 1:modes
   Fourier_sum(i+1,1) = z_in(2*i) + j*z_in(2*i+1); 
   Fourier_sum(modes+i+2) = z_in(2*modes+1+2*i) + j*z_in(2*modes+1+2*i+1); 
end

theta_dot = zeros(2,1);
theta = states(1:2);
theta_tilde = Fourier_sum(1);  
thetadot_tilde = 0;
y_a_ddot = 0;

for k = 1:modes
    theta_tilde = theta_tilde + 2*real(Fourier_sum(k+1)*exp(k*j*2*pi*t));
    thetadot_tilde = thetadot_tilde + 2*real(Fourier_sum(k+1)*(k*j*2*pi)*exp(k*j*2*pi*t));
    y_a_ddot = y_a_ddot + 2*real(Fourier_sum(modes+k+i)*(k*j*2*pi)^2*exp(k*j*2*pi*(t-tau)));                   
end

y_a_ddot = y_a_ddot*(timescale^2); 
thetadot_tilde = thetadot_tilde*timescale;

%dydt = [thetadot; thetaddot; ydot; yddot]
%y = [theta; thetadot; y; ydot]

states_dot = zeros(4,1);
states_dot(1) = states(2);
states_dot(3) = states(4);

LHS_matrix = [m*l^2, 0; m*l*sin(states(1)), M].*timescale^2;
RHS_matrix = [-m*g*l*sin(states(1)) - kappa*states(2)*timescale - m*l*y_a_ddot*sin(states(1)) - k1*(states(1) - theta_tilde) - k2*(states(2)*timescale - thetadot_tilde);
              -C*states(4)*timescale - K*states(3) - m*l*(states(2)^2)*cos(states(1))*(timescale^2) + a*cos(4*pi*t) - m*y_a_ddot];

states_dot([2,4]) = LHS_matrix \ RHS_matrix;
end
