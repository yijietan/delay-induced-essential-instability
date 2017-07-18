function [maxtheta,z_new] = Newton_continuation(z_0, numerical_parameters, physical_parameters,omega,time_array)

runtime = numerical_parameters(1);
points = numerical_parameters(2);
modes = numerical_parameters(3);
stepsize = numerical_parameters(4);
tolerance = numerical_parameters(5);
continuation_stepsize = numerical_parameters(6);
continuation_points = numerical_parameters(7);
direction = numerical_parameters(8);

z_new(:,1) = z_0;
a_normalized = z_0(end);


for i = 1:continuation_points
    
    if i == 1
        z = z_0(1:end-1);
        [J,F] = Jacobian_continuation(z,a_normalized,omega,modes,runtime,points,stepsize,physical_parameters);
        a_norm_delta = a_normalized + 1e-4;              
        [time,states] = delay_system(z,a_norm_delta,omega,modes,runtime,points,physical_parameters);
        z_out = Fourier(time,states,modes,runtime);
        %[theta_recon, y_recon] = test(z_out,modes,time_array,states);
        F_delta = residual_continuation(z,z_out,modes,physical_parameters);
        F_a = (F_delta - F)/1e-4;
%         z_tan = [-J\F_a; -1] / norm([-J\F_a; -1], 2);
        x = -J\F_a;
        alpha_prime = -(1 + x.'*x)^(-1/2);
        z_prime = x*alpha_prime;
        z_tan = [z_prime; alpha_prime];
        z_k = z_0 + continuation_stepsize*z_tan;
        [J_k,F_k,F_a] = Broyden_continuation(J,F,z_0,z_k,modes,runtime,points,stepsize,physical_parameters,omega);
    else
%         z_tan = [-J_k\F_a; direction] / norm([-J_k\F_a; direction], 2);
        x = -J_k\F_a;
        alpha_prime = -(1 + x.'*x)^(-1/2);
        z_prime = x*alpha_prime;
        z_tan = [z_prime; alpha_prime];
        z_k = z_0 + continuation_stepsize*z_tan;
        [J_k,F_k,F_a] = Broyden_continuation(J_k,F_k,z_0,z_k,modes,runtime,points,stepsize,physical_parameters,omega);
    end
    
    F_abs = 1;
    while F_abs > tolerance
        z_2 = (J_k)\(-F_a);
        z_1 = (J_k)\(-F_k);
        g = (z_k(1:end-1) - z_0(1:end-1))'*z_tan(1:end-1) + (z_k(end) - z_0(end))*z_tan(end) - continuation_stepsize;
        
        inc_a_kplus1 = -(g + z_1'*z_tan(1: end-1)) / (z_tan(end) + z_2'*z_tan(1: end-1));
        inc_z_kplus1 = z_1 + z_2*inc_a_kplus1;
        
        z_kplus1 = z_k + [inc_z_kplus1; inc_a_kplus1];
        
        [J_kplus1, F_kplus1, F_a_kplus1] = Broyden_continuation(J_k, F_k, z_k, z_kplus1, modes, runtime, points, stepsize, physical_parameters, omega);
        z_k = z_kplus1;
        J_k = J_kplus1;
        F_k = F_kplus1;
        F_a = F_a_kplus1;
        i
        F_abs = norm(abs(F_kplus1), 1);
    end
    z_0 = z_k;
    z_new(:,i+1) = z_k;
    maxtheta(:,i+1) = max_theta(z_k,modes);
end
end