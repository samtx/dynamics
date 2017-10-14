% HW Assignment 2
% Sam Friedman
% 10/14/2017
% refer to p. 248 in textbook

function [] = hw2()



L1 = 0.241;  % [m]
L2 = 1.206;  % [m]
% Solve the kinematics of the system, assuming a constant theta_dt
theta_dt = 39.5;  % [rad/s], constant angular velocity
dt = 0.001;         % [s],  time step
n = ceil(2*pi/(theta_dt*dt)) ; % total time steps, rounded up to nearest integer
T = 0:dt:(n-1)*dt';
theta = zeros(n+1,1);  % initialize state vectors
phi = zeros(n,1);
phi_dt = zeros(n,1);
phi_dt2 = zeros(n,1);

for i = 1:n
    phi(i) = asin(L1/L2*sin(theta(i)));
    phi_dt(i) = theta_dt*(L1*cos(theta(i)))/(L2*cos(phi(i)));
    phi_dt2(i) = 1/(L2*cos(phi(i)))*(L2*sin(phi(i))*phi_dt(i)^2-L1*sin(theta(i))*theta_dt^2);
    theta(i+1) = dt*theta_dt + theta(i);
end
% remove last theta value
theta(end) = [];
% plot figures
subplot(2,2,1)
plot(T,theta)
xlabel('t [s]');
ylabel('rad');
legend('\theta','Location','best');

subplot(2,2,2)
plot(T,phi)
xlabel('t [s]');
ylabel('rad');
legend('\phi','Location','best');

subplot(2,2,3)
plot(T,phi_dt)
xlabel('t [s]');
ylabel('rad/s');
legend('\dot{\phi}','Location','best','Interpreter','latex');

subplot(2,2,4)
plot(T,phi_dt2);
xlabel('t [s]');
ylabel('rad/s2');
legend('\ddot{\phi}','Location','best','Interpreter','latex');


% A = [
%                Ia,                   0,   -L1*sin(theta),    L1*cos(theta),                0,               0;
%                 0,                 I2g, -0.5*L2*sin(phi), -0.5*L2*cos(phi), -0.5*L2*sin(phi), 0.5*L2*cos(phi);
% -m2*L1*sin(theta), -0.5*m2*L2*sin(phi),    
% ];

end


function dydt = f(t,y,constants)
   g   = constants.g;  % gravitational acceleration
   m1  = constants.m1;
   L1  = constants.L1;
   m2  = constants.m2;
   L2  = constants.L2;
   mp  = constants.mp;
   Ia  = 1/3*m1*L1^2;
   I2g = 1/12*m2*L2^2;
   w1  = m1*g;
   w2  = m2*g;
   dydt = [
       
   ];
   
end