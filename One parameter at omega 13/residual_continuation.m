function F = residual_continuation(z_in,z_out,modes,physical_parameters)
tau = physical_parameters(10);

F = zeros(4*modes+2,1);
for i = 1:2*modes+2
    F(i,1) = z_in(i) - z_out(i);
end

n = 1;
for i = 2*modes+3:2:4*modes+2
    if mod(i,2) ~= 0
        F(i) = z_in(i)*cos(2*pi*n*tau) + z_in(i+1)*sin(2*pi*n*tau) - z_out(i);
        F(i+1) = -z_in(i)*sin(2*pi*n*tau) + z_in(i+1)*cos(2*pi*n*tau) - z_out(i+1);
    end
    n = n + 1;
end