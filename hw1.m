% HW 1

% Constants
lb_to_kg = 0.4536;
mh = 3.4 * lb_to_kg;  %kg
mr = 8.7 * lb_to_kg;  %kg
pAll = [50; 100; 200]; % psi
kzAll = [10; 20; 30]*10^6;  % N/m
czAll = [4; 6; 9]*10^3;   % N.s/m
F = [350; 0];

% mass matrix
M = [mh, 0;
     0 ,mr];

i = 1;
cz = czAll(i);
kz = kzAll(i);

% damping matrix
C = [cz, -cz;
    -cz,2*cz];

% stiffness matrix
K = [kz, -kz;
    -kz,2*kz];

% Find eigenvalues for given kz
[A,e] = eig(K,M);
A = A*-1;
w = sqrt(diag(e));  % get natural frequencies from eigenvalues

% get normalized, uncoupled coeffs from eigenvector matrix A
Cq = A'*C*A;
Mq = A'*M*A;
Kq = A'*K*A;
Fq = A'*F;

fq = diag(Fq);
kq = diag(Kq);

% get damping coeffs
zeta = diag(Cq)./(2.*w);


% solve 2nd order nonhomogeneous ODE for q
q = @(t,y) [y(2);
            fq(i)-2*zeta(i)*w(i)*y(2)-kq(i)*y(1)];
[t, y] = ode45(q,[0,0.016],[0,0]);

plot(t,y(:,1));
xlabel('t');
ylabel('q(t)');