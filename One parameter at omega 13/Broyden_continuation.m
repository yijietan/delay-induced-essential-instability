function [J_kplus1, F_kplus1, F_a] = Broyden_continuation(J_k, F_k, z_k, z_kplus1, modes, runtime, points, continuation_stepsize, physical_parameters, omega)

a_normalized = z_kplus1(end);

z_in = z_kplus1(1:end-1);
[time,states] = delay_system(z_in,a_normalized,omega,modes,runtime,points,physical_parameters);
z_out = Fourier(time,states,modes,runtime);

F_kplus1 = residual_continuation(z_in,z_out,modes,physical_parameters);
J_kplus1 = J_k + ((F_kplus1-F_k) - J_k*(z_kplus1(1:end-1) - z_k(1:end-1)))*(z_kplus1(1:end-1) - z_k(1:end-1))' / ((z_kplus1(1:end-1) - z_k(1:end-1))'*(z_kplus1(1:end-1) - z_k(1:end-1)));

% Obtaining partial derivative F_a
a_normalized_delta = a_normalized + continuation_stepsize;
[time,states] = delay_system(z_in,a_normalized_delta,omega,modes,runtime,points,physical_parameters);
z_out = Fourier(time,states,modes,runtime);
F_delta_a = residual_continuation(z_in,z_out,modes,physical_parameters);
F_a = (F_delta_a - F_kplus1) / continuation_stepsize;
end