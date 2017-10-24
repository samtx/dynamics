% HW Assignment 2
% Sam Friedman, Benson Isaac, Mohamed Mohamed, Alexis Trevino
% 10/24/2017

function [] = hw2()

% Constants for problem
L1 = 0.241;         % [m]      length of arm AB
L2 = 1.206;         % [m]      length of arm BC
m1 = 908;           % [kg]     mass of arm AB
m2 = 245;           % [kg]     mass of arm BC
mp = 590;           % [kg]     mass of piston
g = 9.81;           % [m/s^2]  gravitational constant
w1 = m1*g;          % [N]      weight of arm AB
w2 = m2*g;          % [N]      weight of arm BC
I2g = m2*L2^2/12;   % [kg.m^2] moment of inertia

% if it doesn't exist, create data directory to save results in.
if ~isdir('data')
    mkdir('data');
end

% Solve the kinematics of the system, assuming a constant theta_dt
theta_dt_list = [39.5, 79.0];      % angular velocity in rad/s
rpm_list = {'377rpm','754rpm'};    % string of RPM values for filenames for saving data

% Loop over each given constant value of theta_dt
for j = 1:length(theta_dt_list)
    theta_dt = theta_dt_list(j);   % [rad/s], constant angular velocity
    dt = 0.001;                    % [s],  time step
    n = ceil(2*pi/(theta_dt*dt)) ; % total time steps for one revolution, rounded up to nearest integer
    t = [0:dt:(n-1)*dt]'; %#ok<NBRAK>
    theta = t*theta_dt;  % vector of theta values
    
    % Initialize state vectors
    phi = zeros(n,1);
    phi_dt = zeros(n,1);
    phi_dt2 = zeros(n,1);
    xp = zeros(n,1);
    xp_dt = zeros(n,1);
    xp_dt2 = zeros(n,1);
    
    % Compute values for phi, phi_dt, phi_dt2, xp, xp_dt, xp_dt2 for each value of theta
    for i = 1:n
        phi(i) = asin(L1/L2*sin(theta(i)));
        phi_dt(i) = theta_dt*(L1*cos(theta(i)))/(L2*cos(phi(i)));
        phi_dt2(i) = 1/(L2*cos(phi(i)))*(L2*sin(phi(i))*phi_dt(i)^2-L1*sin(theta(i))*theta_dt^2);
        xp(i) = L1*cos(theta(i))+L2*cos(phi(i));
        xp_dt(i) = -L1*sin(theta(i))*theta_dt-L2*sin(phi(i))*phi_dt(i);
        xp_dt2(i) = -L1*cos(theta(i))*theta_dt^2-L2*(cos(phi(i))*phi_dt(i)^2+sin(phi(i))*phi_dt2(i));
    end
    
    % Initialize moment and reaction force vectors
    M = zeros(n,1);   % [N.m]    applied moment
    Bx = zeros(n,1);  % [N]      reaction force in x-direction at location B
    By = zeros(n,1);  % [N]      reaction force in y-direction at location B
    Cx = zeros(n,1);  % [N]      reaction force in x-direction at location C
    Cy = zeros(n,1);  % [N]      reaction force in y-direction at location C
    
    % For the each value of theta, compute values for the applied moment M and reaction forces Bx, By, Cx, Cy
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
        
        x = A\b;       % solve matrix equation
        M(i)  = x(1);  % store data values
        Bx(i) = x(2);
        By(i) = x(3);
        Cx(i) = x(4);
        Cy(i) = x(5);
    end
    
    % Compute Ax, Ay reaction forces
    Ax = -Bx;
    Ay = w1 + By;
    
    % Compute net reaction forces
    Anet = sqrt(Ax.^2+Ay.^2); %#ok<NASGU>
    Bnet = sqrt(Bx.^2+By.^2); %#ok<NASGU>
    Cnet = sqrt(Cx.^2+Cy.^2); %#ok<NASGU>
    
    % Save data to .dat files. Plots are made in LaTeX document report using Pgfplots.
    % Create one file for each variable. Include the time steps in vector t.
    % Store files in data folder
    vars = {'theta','phi','phi_dt','phi_dt2','xp','xp_dt','xp_dt2','M','Bx','By','Cx','Cy','Ax','Ay','Anet','Bnet','Cnet'};
    for i = 1:length(vars)  % iterate through variables to save data
        tbl = array2table([t, eval(vars{i})]);
        writetable(tbl,['data/',rpm_list{j},'_',vars{i},'.dat'],'WriteVariableNames',false,'Delimiter',' ');
    end
    
end

end