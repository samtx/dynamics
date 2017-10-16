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
t = [0:dt:(n-1)*dt]'; %#ok<NBRAK>
theta = t*theta_dt;  % theta vector
% initialize state vectors
phi = zeros(n,1);
phi_dt = zeros(n,1);
phi_dt2 = zeros(n,1);
xp = zeros(n,1);
xp_dt = zeros(n,1);
xp_dt2 = zeros(n,1);

for i = 1:n
    phi(i) = asin(L1/L2*sin(theta(i)));
    phi_dt(i) = theta_dt*(L1*cos(theta(i)))/(L2*cos(phi(i)));
    phi_dt2(i) = 1/(L2*cos(phi(i)))*(L2*sin(phi(i))*phi_dt(i)^2-L1*sin(theta(i))*theta_dt^2);
    xp(i) = L1*cos(theta(i))+L2*cos(phi(i));
    xp_dt(i) = -L1*sin(theta(i))*theta_dt-L2*sin(phi(i))*phi_dt(i);
    xp_dt2(i) = -L1*cos(theta(i))*theta_dt^2-L2*(cos(phi(i))*phi_dt(i)^2+sin(phi(i))*phi_dt2(i));
end

m1 = 908;    % [kg]
m2 = 245;    % [kg]
mp = 590;
g = 9.81;    % [m/s^2]  gravitational constant
w1 = m1*g; % [N]
w2 = m2*g; % [N]
I2g = m2*L2^2/12; % [m.o.i.]
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

% Compute Ax, Ay reaction forces
Ax = -Bx; %#ok<NASGU>
Ay = w1 + By; %#ok<NASGU>

% Save data to .dat files
vars = {'theta','phi','phi_dt','phi_dt2','xp','xp_dt','xp_dt2','M','Bx','By','Cx','Cy','Ax','Ay'};
for i = 1:length(vars)  % iterate through variables to save data
	tbl = array2table([t, eval(vars{i})]);
    writetable(tbl,['data/',vars{i},'.dat'],'WriteVariableNames',false,'Delimiter',' ');
end

% skip plotting for now...
return

% plot figures
figure(1)
subplot(2,2,1)
plot(t,theta)
xlabel('t [s]');
ylabel('rad');
legend('\theta','Location','best');

subplot(2,2,2)
plot(t,phi)
xlabel('t [s]');
ylabel('rad');
legend('\phi','Location','best');

subplot(2,2,3)
plot(t,phi_dt)
xlabel('t [s]');
ylabel('rad/s');
legend('$\dot{\phi}$','Location','best','Interpreter','latex');

subplot(2,2,4)
plot(t,phi_dt2);
xlabel('t [s]');
ylabel('rad/s2');
legend('$$\ddot{\phi}$$','Location','best','Interpreter','latex');

figure(2)
subplot(2,2,1)
plot(t,xp)
xlabel('t [s]');
ylabel('m');
legend('xp','Location','best');

subplot(2,2,2)
plot(t,xp_dt)
xlabel('t [s]');
ylabel('m/s');
legend('xp_dt','Location','best');

subplot(2,2,3)
plot(t,xp_dt2)
xlabel('t [s]');
ylabel('m/s^2');
legend('xp_dt2','Location','best','Interpreter','latex');

% plot figures
figure(3)
subplot(2,3,1)
plot(T,M)
xlabel('t [s]');
ylabel('M');
legend('M(t)','Location','best');

hold on
subplot(2,3,2)
plot(T,Bx)
xlabel('t [s]');
ylabel('B_x [N]');
legend('\phi','Location','best');

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