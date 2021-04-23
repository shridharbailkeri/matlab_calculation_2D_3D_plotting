% Program to calculate heat transfer in smooth pipe (one phase)
% Refrigerant; R407C
clear all
% Material properties
Pr = 3.115;
Rho = 1/(0.84*10^-3); % specific volume
% kg / m^3
% dynamic viscosity eta / density = kinematic viscosity
eta = 0.213*10^-3;    % Ns/m^2
lambda = 0.094;       % W / mK

% geometric properties pipe diameter and length

L=1:1:10;             % Length of the profile
L_max = numel(L);
d_i = 1:0.1:10;         % starting from 1 with step size 1 vector is created
d_i_max = numel(d_i);
[d_i, L] = meshgrid(d_i, L);
d_i = d_i .* 10^-3;    % element wise 

A_sa = pi.*d_i.^2./4;   % flow sufrace area
A_ca = pi.*d_i.*L;      % Circumference area (heat transfer)

% Process quantities 
M = 0.03;               % mass flow (kg/s)
% M = [kg/s] = Asa*v*rho = [m^2]*[m/s]*[ kg / m^3]
% v = M / Asa*rho
v = M ./ (A_sa.*Rho);    % velocity
Re = Rho.*d_i.*v./eta;   % Reynolds number
DeltaT = 2;           % 2K temperature difference

% calculation of Heat Transfer loop


for i= 1:L_max
for j=1:d_i_max

switch true
    
    case Re(i,j) > 10000  
    Nu(i,j) = 0.024*Re(i,j)^0.8*Pr^(1/3);
    
    case 10000 > Re(i,j) > 2300
    Xi(i,j) = (1.8*log10(Re(i,j))-1.5)^-2; 
    Nu(i,j) = (Xi(i,j)/8*Re(i,j)*Pr)/(1+12.7*(Xi(i,j)/8)^(0.5)*(Pr^(2/3)-1))*(1+(d_i(i,j)/L)^(2/3));
    
    case Re(i,j) < 2300
    X(i,j) = L/(d_i(i,j)*Re(i,j)*Pr);
    Nu_O(i,j) = 3.657/(tanh(2.264*X(i,j)^(1/3)+1.7*X(i,j)^(2/3)))+0.0499/X(i,j)*tanh(X(i,j));
    Nu(i,j) = Nu_O(i,j) * 1 /(tanh(2.43*Pr^(1/6)*X(i,j)^6));
    
end

alpha(i,j) = Nu(i,j) * lambda/d_i(i,j)     % W/m^K
Q(i,j) = alpha(i,j) * A_ca(i,j) * DeltaT   % W

end
end

mesh(d_i, L, Q)
xlabel('Diameter (m)')
ylabel('Length of the profile in meters')
zlabel('Heat Transfer in W')
title('Heat Transfer Depending on Diameters')
    
    
    
    
        
                       
                       









