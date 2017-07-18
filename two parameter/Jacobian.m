function [J, F] = Jacobian(z,a_normalized,omega,modes,runtime,points,stepsize,physical_parameters, initial_conditions, time_array)
J = zeros(4*modes+3,4*modes+3);  

[time,states] = delay_system(z,a_normalized,omega,modes,runtime,points,physical_parameters);
% [time,states] = ode45(@(t, states)delay_system_test(t, states, physical_parameters, a_normalized, omega, z, modes), time_array, initial_conditions);
% time = time.';
% states = states.';
z_out = Fourier(time,states,modes,runtime);
F = residual(z,z_out,modes,physical_parameters);
z_delta = z;

for i = 1:(length(z) + 1)
    if i < length(z) + 1
        z_delta(i) = z(i) + stepsize;
        a_normalized = a_normalized;
    else
        a_normalized = a_normalized + 1e-5;
    end
    [time,states] = delay_system(z_delta,a_normalized,omega,modes,runtime,points,physical_parameters);
%     [time,states] = ode45(@(t, states)delay_system_test(t, states, physical_parameters, a_normalized, omega, z_delta, modes), time_array, initial_conditions);
%     time = time.';
%     states = states.';
    z_out = Fourier(time,states,modes,runtime);
    F_delta = residual(z_delta,z_out,modes,physical_parameters);  
    J(:,i) = (F_delta - F)/(1e-5);
    z_delta = z;
end

end