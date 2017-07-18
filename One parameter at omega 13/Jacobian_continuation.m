function [J, F] = Jacobian_continuation(z,a_normalized,omega,modes,runtime,points,stepsize,physical_parameters)
%========Caclulating Jacobian for the iteration
J= zeros(4*modes+2,4*modes+2);      %====N: number of modes

[time,states] = delay_system(z,a_normalized,omega,modes,runtime,points,physical_parameters);
z_out = Fourier(time,states,modes,runtime);
F = residual_continuation(z,z_out,modes,physical_parameters);
z_delta = z;

for i = 1:length(z)
    z_delta(i) = z(i) + stepsize;
    [time,states] = delay_system(z_delta,a_normalized,omega,modes,runtime,points,physical_parameters);
    z_out = Fourier(time,states,modes,runtime);
    F_delta = residual_continuation(z_delta,z_out,modes,physical_parameters);  
    J(:,i) = (F_delta - F)/stepsize;
    z_delta = z;
end

end