function [J_kplus1, F_kplus1, F_omega] = Broyden(J_k, F_k, z_k, z_kplus1, modes, runtime, points, stepsize, physical_parameters, initial_conditions, time_array)

omega = z_kplus1(end);
a_normalized = z_kplus1(end-1);

z_in = z_kplus1(1:end-2);
[time,states] = delay_system(z_in,a_normalized,omega,modes,runtime,points,physical_parameters);
% [time,states] = ode45(@(t, states)delay_system_test(t, states, physical_parameters, a_normalized, omega, z_in, modes), time_array, initial_conditions);
% time = time.';
% states = states.';
z_out = Fourier(time,states,modes,runtime);

F_kplus1 = residual(z_in,z_out,modes,physical_parameters);
J_kplus1 = J_k + ((F_kplus1-F_k) - J_k*(z_kplus1(1:end-1) - z_k(1:end-1)))*(z_kplus1(1:end-1) - z_k(1:end-1))' / ((z_kplus1(1:end-1) - z_k(1:end-1))'*(z_kplus1(1:end-1) - z_k(1:end-1)));

% Obtaining partial derivative F_omega
omega_delta = omega + stepsize;
[time,states] = delay_system(z_in,a_normalized,omega_delta,modes,runtime,points,physical_parameters);
% [time,states] = ode45(@(t, states)delay_system_test(t, states, physical_parameters, a_normalized, omega_delta, z_in, modes), time_array, initial_conditions);
% time = time.';
% states = states.';
z_out = Fourier(time,states,modes,runtime);
F_delta_omega = residual(z_in,z_out,modes,physical_parameters);
F_omega = (F_delta_omega - F_kplus1) / stepsize;
end