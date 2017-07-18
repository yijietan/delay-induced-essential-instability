function z_new = Newton_continuation(z_0, numerical_parameters, physical_parameters, initial_conditions, time_array)

runtime = numerical_parameters(1);
points = numerical_parameters(2);
modes = numerical_parameters(3);
stepsize = numerical_parameters(4);
tolerance = numerical_parameters(5);
continuation_stepsize = numerical_parameters(6);
continuation_points = numerical_parameters(7);
direction = numerical_parameters(8);

z_new(:,1) = z_0;
omega = z_0(end);
a_normalized = z_0(end-1);


for i = 1:continuation_points
    
    if i == 1
        z = z_0(1:end-2);
        [J,F] = Jacobian(z,a_normalized,omega,modes,runtime,points,stepsize,physical_parameters, initial_conditions, time_array);
        
        % Find F_omega
        omega_delta = omega + stepsize;
        [time,states] = delay_system(z,a_normalized,omega_delta,modes,runtime,points,physical_parameters);
%         [time,states] = ode45(@(t, states)delay_system_test(t, states, physical_parameters, a_normalized, omega_delta, z, modes), time_array, initial_conditions);
%         time = time.';
%         states = states.';
        z_out = Fourier(time,states,modes,runtime);
        F_delta = residual(z,z_out,modes,physical_parameters);
        F_omega = (F_delta - F)/stepsize;
        
        % Calculate tangent vector
        z_tan = [-J\F_omega; 1] / norm([-J\F_omega; 1], 2);
        z_k = z_0 + continuation_stepsize*z_tan;
        [J_k,F_k,F_omega] = Broyden(J,F,z_0,z_k,modes,runtime,points,stepsize,physical_parameters,initial_conditions,time_array);
    else
        z_tan = [-J_k\F_omega; direction] / norm([-J_k\F_omega; direction], 2);
        z_k = z_0 + continuation_stepsize*z_tan;
        [J_k,F_k,F_omega] = Broyden(J_k,F_k,z_0,z_k,modes,runtime,points,stepsize,physical_parameters,initial_conditions,time_array);
    end
    
    F_abs = 1;
    while F_abs > tolerance
        % Carry out pseudo-arclength algorithm
        z_2 = (J_k)\(-F_omega);
        z_1 = (J_k)\(-F_k);
        g = (z_k(1:end-1) - z_0(1:end-1))'*z_tan(1:end-1) + (z_k(end) - z_0(end))*z_tan(end) - continuation_stepsize;
        
        inc_omega_kplus1 = -(g + z_1'*z_tan(1: end-1)) / (z_tan(end) + z_2'*z_tan(1: end-1));
        inc_z_kplus1 = z_1 + z_2*inc_omega_kplus1;
        
        z_kplus1 = z_k + [inc_z_kplus1; inc_omega_kplus1];
        
        [J_kplus1, F_kplus1, F_omega_kplus1] = Broyden(J_k, F_k, z_k, z_kplus1, modes, runtime, points, stepsize, physical_parameters, initial_conditions, time_array);
        z_k = z_kplus1;
        J_k = J_kplus1;
        F_k = F_kplus1;
        F_omega = F_omega_kplus1;
        i
        F_abs = norm(abs(F_kplus1), 1);
    end
    z_0 = z_k;
    z_new(:,i+1) = z_k;
end
end