function maxtheta = max_theta(z,modes)
theta = z(1);
t = linspace(0,1,1000);
dt = t(end) - t(end-1);
for i = 1:modes
    theta_complex = z(2*i) + j*z(2*i+1);
    theta = theta + 2*real(theta_complex*exp(i*j*2*pi.*(t)));
end
maxtheta = max(theta);
end