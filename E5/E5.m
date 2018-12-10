%%
tmp = matlab.desktop.editor.getActive; %% cd to current path
cd(fileparts(tmp.Filename));
addpath('../H2/m_scripts')

%% task 1
N = 201;
r_max = 5;
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

%% task 3: minimize E instead of epsilon
clc;
clear; clf;

r_max = 8;
N = (40*r_max)+1;
r_min = 0;
dr = (r_max-r_min)/(N-1);
rr = linspace(r_min,r_max,N)';
rrLarge =  rr * ones(1,N);

% build matrices
A = ( diag(ones(N-1,1),1) + diag(ones(N-1,1),-1) + -2*eye(N) );
A = A/dr^2;
A(1,1) = 1;
A(1,2) = 0;
A(end,end) = 1;
A(end, end-1) = 0;

for Z = [2 27/16];
    H0 = -0.5*A-diag(Z./rr); % 2 = He
    H0(1,1) = 0; % makes epsilon * f_1  = 0;
    H0(end,end) = 0; % makes epsilon * f_N  = 0;
    V = zeros(size(rr));
    dE = 1;E_old = 1; n = 0;
    
    while abs(dE) > 1e-4
        H = H0 + diag([0; V(2:end-1); 0]);
        [F,D] = eig(H);
        epsilon = diag(D);
        ns_vec = F(:,:).^2./ (4*pi.*rrLarge.^2)/dr ; % all eigenvectors
        ns_vec(1,:) = deal(0);
        
        rhs_plus_BC = [zeros(size(ns_vec(1,:)));
            -4*pi*rrLarge(2:end-1,:) .* ns_vec(2:end-1,:);
            ones(size(ns_vec(1,:)))];
        U = A\rhs_plus_BC;
        V = U./rrLarge;
        V(1,:) = deal(1);
        
        E0s = 2*epsilon' - 4*pi*trapz(rr, rrLarge.^2.*V.*ns_vec, 1);
        
        [E0,i] = min(E0s);
        V = V(:,i);
        
        dE = E0-E_old;
        
        E_old = E0;
        n = n+1;
        fprintf('E = %.4f\n', E0)
        
    end
    fprintf('n = %d\t', n);
    fprintf('E = %.4f\n', E0)
    n_s = ns_vec(:,i);
    if Z ==2
        plot(rr, 4*pi*rr.^2.*n_s, '-r','linewidth',2); hold on;
        plot(rr, 4*rr.^2.*Z^3.*exp(-2*Z*rr), ':r','linewidth',2); hold on;
    else
        plot(rr, 4*pi*rr.^2.*n_s, '-k','linewidth',2); hold on;
        plot(rr, 4*rr.^2.*Z^3.*exp(-2*Z*rr), ':k','linewidth',2); hold on;
    end
end

legend('Z = 2, Hartree', 'Z = 2, central-field approx', 'Z = 27/16, Hartree', 'Z = 27/16, central-field approx')

xlabel('r [a_0]');
ylabel('\rho(r)');




