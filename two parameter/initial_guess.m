function z_initialguess = initial_guess(runtime, physical_parameters, omega, a_normalized, modes, initial_conditions, time_array)
[~,states] = ode45(@(t, states)nodelay_system(t, states, physical_parameters, a_normalized, omega), time_array, initial_conditions);
states = states.';
z_initialguess = Fourier(time_array, states, modes, runtime);
end