function [fourier_coeff]=Fourier(time_array,states,modes,runtime)
theta = states(1,:);
y = states(3,:);

period = 1;
dt = time_array(end) - time_array(end-1);
n1 = 1 + floor((runtime-period)/dt); 
n2 = length(time_array);
y_oneperiod = y(:,n1:n2); 
theta_oneperiod = theta(:,n1:n2);
t = time_array(:,n1:n2);

thetaFFT = fft(theta_oneperiod/length(theta_oneperiod));
yFFT = fft(y_oneperiod/length(y_oneperiod));


fourier_coeff=zeros(4*modes+2,1);
for i=1:modes+1
    if i == 1
      fourier_coeff(i) = real(thetaFFT(i));
      fourier_coeff(i+2*modes+1) = real(yFFT(i));
    else
      fourier_coeff(2*(i-1),1)= real(thetaFFT(i));
      fourier_coeff(2*(i-1)+1,1)= imag(thetaFFT(i));
       
      fourier_coeff(2*(i-1) + 2*modes + 1, 1)= real(yFFT(i));
      fourier_coeff(2*(i-1) + 2*modes + 2, 1)= imag(yFFT(i));
    end
end
end