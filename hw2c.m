% HW Assignment 2
% Sam Friedman
% 10/14/2017
% refer to p. 248 in textbook

function [] = hw2c()

% get state vector data from previous function
load('hw2_data.mat');

% L1 = 0.241;  % [m]
% L2 = 1.206;  % [m]
m1 = 908;    % [kg]
m2 = 245;    % [kg]
mp = 590;
g = 9.81;    % [m/s^2]  gravitational constant
w1 = m1*g; % [N]
w2 = m2*g; % [N]
I2g = m2*L2^2/12; % [m.o.i.]
% Solve the kinematics of the system, assuming a constant theta_dt
% theta_dt = 39.5;  % [rad/s], constant angular velocity
% dt = 0.001;         % [s],  time step
% n = ceil(2*pi/(theta_dt*dt)) ; % total time steps, rounded up to nearest integer
% T = 0:dt:(n-1)*dt';
% theta = zeros(n+1,1);  % initialize state vectors
M = zeros(n,1); % applied moment
Bx = zeros(n,1);  % reaction forces
By = zeros(n,1);
Cx = zeros(n,1);
Cy = zeros(n,1);

for i = 1:n
    A = [
    1,    L1*sin(theta(i)),   -L1*cos(theta(i)),                   0,                  0;
    0, -0.5*L2*sin(phi(i)), -0.5*L2*cos(phi(i)), -0.5*L2*sin(phi(i)), 0.5*L2*cos(phi(i));
    0,                  -1,                   0,                   1,                  0;    
    0,                   0,                  -1,                   0,                 -1;
    0,                   0,                   0,                  -1,                  0;
    ];

    b = [
    w1*0.5*L1*cos(theta(i));
    -I2g*phi_dt2(i);
    m2*L1*theta_dt^2*cos(theta(i))+m2*0.5*L2*cos(phi(i))*phi_dt(i)^2+m2*0.5*L2*sin(phi(i))*phi_dt2(i);
    m1*0.5*L1*theta_dt^2*sin(theta(i))-w2;
    mp*L1*theta_dt^2*cos(theta(i))+mp*L2*phi_dt(i)^2*cos(phi(i))+mp*L2*sin(phi(i))*phi_dt2(i);
    ];

    x = A\b;  % solve matrix equation
    M(i) = x(1); % store data values
    Bx(i) = x(2);
    By(i) = x(3);
    Cx(i) = x(4);
    Cy(i) = x(5);
     
end

save('hw2_data2.mat') % save data to .mat file

% plot figures
figure(1)
subplot(2,3,1)
plot(T,M)
xlabel('t [s]');
ylabel('M');
legend('M(t)','Location','best');

subplot(2,3,2)
plot(T,Bx)
xlabel('t [s]');
ylabel('B_x [N]');
% legend('\phi','Location','best');

subplot(2,3,3)
plot(T,By)
xlabel('t [s]');
ylabel('B_y');
% legend('$\dot{\phi}$','Location','best','Interpreter','latex');

subplot(2,3,4)
plot(T,Cx);
xlabel('t [s]');
ylabel('C_x');
% legend('$$\ddot{\phi}$$','Location','best','Interpreter','latex');

subplot(2,3,5)
plot(T,Cy)
xlabel('t [s]');
ylabel('C_y');
% legend('xp','Location','best');

end