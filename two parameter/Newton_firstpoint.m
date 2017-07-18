function z_firstpoint = Newton_firstpoint(z_initialguess, a_normalized, omega, numerical_parameters, physical_parameters, initial_conditions, time_array)
runtime = numerical_parameters(1);
points = numerical_parameters(2);
modes = numerical_parameters(3);
stepsize = numerical_parameters(4);
tolerance = numerical_parameters(5);

% Calculate Jacobian for first step
J = zeros(4*modes+3,4*modes+3);
[time,states] = delay_system(z_initialguess,a_normalized,omega,modes,runtime,points,physical_parameters);
% [time,states] = ode45(@(t, states)delay_system_test(t, states, physical_parameters, a_normalized, omega, z_initialguess, modes), time_array, initial_conditions);
% time = time.';
% states = states.';
z = Fourier(time,states,modes,runtime);
F_l = residual(z_initialguess,z,modes,physical_parameters);
z_initialguess_delta = z_initialguess;

for i = 1:4*modes+3
    if i ~= 4*modes+3
        z_initialguess_delta(i) = z_initialguess(i) + stepsize;
        a_normalized = a_normalized;
    else
        a_normalized = a_normalized + 1e-5;
    end
    [time,states] = delay_system(z_initialguess_delta,a_normalized,omega,modes,runtime,points,physical_parameters);
%     [time,states] = ode45(@(t, states)delay_system_test(t, states, physical_parameters, a_normalized, omega, z_initialguess, modes), time_array, initial_conditions);
%     time = time.';
%     states = states.';
    z_delta = Fourier(time,states,modes,runtime);
    F_delta = residual(z_initialguess_delta,z_delta,modes,physical_parameters);  
    J(:,i) = (F_delta - F_l)/stepsize;
    z_initialguess_delta = z_initialguess;
end

z_l = [z_initialguess;a_normalized];

F_abs = 1;

while F_abs > tolerance
    z_lplus1 = z_l - J\F_l;
    z_in = z_lplus1(1:end-1);
    a_normalized = z_lplus1(end);
    [time,states] = delay_system(z_in,a_normalized,omega,modes,runtime,points,physical_parameters);
%     [time,states] = ode45(@(t, states)delay_system_test(t, states, physical_parameters, a_normalized, omega, z_in, modes), time_array, initial_conditions);
%     time = time.';
%     states = states.';
    z_out = Fourier(time,states,modes,runtime);
    F_lplus1 = residual(z_in,z_out,modes,physical_parameters);
    F_abs = norm(abs(F_lplus1), 1)
    
    % Broyden rank 1 update
    J = J + ((F_lplus1 - F_l) - J*(z_lplus1 - z_l))*(z_lplus1 - z_l)' / ((z_lplus1 - z_l)'*(z_lplus1 - z_l));
    F_l = F_lplus1;
    z_l = z_lplus1;
    z_firstpoint = z_l;
end
end