function F = residual(z_in,z_out,modes,physical_parameters)
tau = physical_parameters(10);
r = physical_parameters(11);

F = zeros(4*modes+3,1);
for i = 1:2*modes+2
    F(i,1) = z_in(i) - z_out(i);
end

k = 1;
for i = 2*modes+3:2:4*modes+2
    if mod(i,2) ~= 0
        F(i) = z_in(i)*cos(2*pi*k*tau) + z_in(i+1)*sin(2*pi*k*tau) - z_out(i);
        F(i+1) = -z_in(i)*sin(2*pi*k*tau) + z_in(i+1)*cos(2*pi*k*tau) - z_out(i+1);
    end
    k = k + 1;
end

% Small amplitude condition
theta = z_in(1);
t = linspace(0,1,1000);
dt = t(end) - t(end-1);
for i = 1:modes
    theta_complex = z_in(2*i) + j*z_in(2*i+1);
    theta = theta + 2*real(theta_complex*exp(i*j*2*pi.*(t)));
end

F(end) = trapz(theta.^2)*dt - r^2;
end