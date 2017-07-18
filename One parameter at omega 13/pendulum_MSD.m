function dydt = pendulum_MSD(t,y,parameter,omega)  

m = 0.27; 
l = 0.1955; 
kappa = 7.5e-3;
k1 = 4; 
k2 = 4; 
M = 0.2; 
C = 0.1; 
K = 200; 
g = 9.81; 
a = parameter*sqrt(omega^2*C^2 + K^2);
timescale = omega / (4*pi); %to scale outputs to have a period of 1

%dydt = [thetadot; thetaddot; ydot; yddot]
%y = [theta; thetadot; y; ydot]

dydt=zeros(4,1);
dydt(1) = y(2);
dydt(3) = y(4);

LHS_matrix = [m*l^2, m*l*sin(y(1)); m*l*sin(y(1)), M+m].*timescale^2;
RHS_matrix = [-kappa*y(2)*timescale - m*g*l*sin(y(1));
              -C*y(4)*timescale - K*y(3) - m*l*(y(2)^2)*cos(y(1))*(timescale^2) + a*cos(4*pi*t)];
          
dydt([2,4])=LHS_matrix \ RHS_matrix;
end