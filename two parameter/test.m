function theta_recon = test(z,modes,time_array,states)
theta = states(1,:);
y = states(3,:);

theta_recon = z(1);
for i = 1:modes
    theta_complex = z(2*i) + j*z(2*i+1);
    theta_recon = theta_recon + 2*real(theta_complex*exp(i*j*2*pi.*(time_array)));
    
    y_complex = z(2*i + 2*modes+1) + j*z(2*i+1 + 2*modes+1);
    y_recon = y_complex + 2*real(y_complex*exp(i*j*2*pi.*(time_array)));
end

figure;
plot(time_array, theta);
hold all;
plot(time_array, theta_recon);
legend('theta, theta_recon');

figure;
plot(time_array, y);
hold all;
plot(time_array, y_recon);
legend('y, y_recon');

end