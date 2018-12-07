%% task 1
N = 201;
r_max = 10;
r_min = 0;

dr = (r_max-r_min)/(N-1);
rr = linspace(r_min,r_max,N)';

A = ( diag(ones(N-1,1),1) + diag(ones(N-1,1),-1) + -2*eye(N) );
A = A/dr^2;
A(1,1) = 1;
A(1,2) = 0;
A(end,end) = 1;
A(end, end-1) = 0;

n_s = 1/pi * exp(-2*rr); % units of Bohr radii
rhs_plus_BC = -4*pi*rr .* n_s;
rhs_plus_BC(1) = 0;
rhs_plus_BC(end) = 1;

U_test = A\rhs_plus_BC;

V_Hartree = 1./rr - (1+ 1./rr).*exp(-2*rr);
U_Hartree = rr.* V_Hartree;

clf;
plot(rr, U_test , '.-',rr, U_Hartree, ':k')
xlabel('r [a_0]');
ylabel('U');

%% task 2

H0 = -0.5*A-diag(1./rr); % the matrix
H0(1,1) = 0; % makes epsilon * f_1  = 0;
H0(end,end) = 0; % makes epsilon * f_N  = 0;
[F,D] = eig(H0);
clf;


[emin,i] = min(min(D));
E0 = D(i,i)
plot(rr, F(:,i)./(sqrt(4*pi).*rr)/sqrt(dr),'r', rr, sqrt(n_s), ':k','linewidth',2)
xlabel('r [a_0]');
ylabel('n(r)');
legend('numerical', 'theoretical')

%% task 3
clc;
clear; clf;
N = 201;
r_max = 5;
r_min = 0;
dr = (r_max-r_min)/(N-1);
rr = linspace(r_min,r_max,N)';

% build matrices
A = ( diag(ones(N-1,1),1) + diag(ones(N-1,1),-1) + -2*eye(N) );
A = A/dr^2;
A(1,1) = 1;
A(1,2) = 0;
A(end,end) = 1;
A(end, end-1) = 0;

Z = 2;
H0 = -0.5*A-diag(Z./rr); % 2 = He
H0(1,1) = 0; % makes epsilon * f_1  = 0;
H0(end,end) = 0; % makes epsilon * f_N  = 0;
V = zeros(size(rr));

dE = 1;E_old = 1; n = 0;
while abs(dE) > 1e-4
    H = H0 + diag([0; V(2:end-1); 0]);
    [F,D] = eig(H);
    [epsilon,i] = min(min(D));
    n_s = F(:,i).^2./(4*pi.*rr.^2)/dr;
    n_s(1) = 0;
    
    rhs_plus_BC = [0; n_s(2:end-1); 1];
    U = A\rhs_plus_BC;
    V = U./rr;
    V(1) = 1;
    
    E0 = 2*epsilon - 4*pi*trapz(rr, rr.^2.*V.*n_s); % this might be wrong. normalization?
    % test normalization. 4*pi*trapz(rr, rr.^2.*V.*n_s)
    dE = E0-E_old;
    
    E_old = E0;
    n = n+1;
    
    fprintf('E = %.4f\n', E0)
end
fprintf('n = %d\n\n', n);
plot(rr, 4*pi*rr.^2.*n_s, '-r'); hold on;
plot(rr, 4*rr.^2.*Z^3.*exp(-2*Z*rr), ':r'); hold on;

Z = 27/16;
H0 = -0.5*A-diag(Z./rr); % 2 = He
H0(1,1) = 0; % makes epsilon * f_1  = 0;
H0(end,end) = 0; % makes epsilon * f_N  = 0;
%V = zeros(size(rr));
dE = 1;E_old = 1; n = 0;
while abs(dE) > 1e-4
    H = H0 + diag([0; V(2:end-1); 0]);
    [F,D] = eig(H);
    [epsilon,i] = min(min(D));
    n_s = F(:,i).^2./(4*pi.*rr.^2)/dr;
    n_s(1) = 0;
    
    rhs_plus_BC = [0; n_s(2:end-1); 1];
    U = A\rhs_plus_BC;
    V = U./rr;
    V(1) = 1;
    
    E0 = 2*epsilon - 4*pi*trapz(rr, rr.^2.*V.*n_s);
    dE = E0-E_old;
    
    E_old = E0;
    n = n+1;
    
    fprintf('E = %.4f\n', E0)
end
fprintf('n = %d\n\n', n);
plot(rr, 4*pi*rr.^2.*n_s, '-k'); hold on;
plot(rr, 4*rr.^2.*Z^3.*exp(-2*Z*rr), ':k'); hold on;
legend('Z = 2, Hartree', 'Z = 2, central-field approx', 'Z = 27/16, Hartree', 'Z = 27/16, central-field approx')

xlabel('r [a_0]');
ylabel('\rho(r)');
