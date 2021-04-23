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

d_i = 1:0.1:10          % starting from 1 with step size 1 vector is created
d_i = d_i .* 10^-3;   % element wise 
L = 2;                % meters 
A_sa = pi.*d_i.^2./4;    % flow sufrace area
A_ca = pi.*d_i.*L;      % Circumference area (heat transfer)

% Process quantities 
M = 0.03;             % mass flow (kg/s)
% M = [kg/s] = Asa*v*rho = [m^2]*[m/s]*[ kg / m^3]
% v = M / Asa*rho
v = M ./ (A_sa.*Rho);   % velocity
Re = Rho.*d_i.*v./eta;   % Reynolds number
DeltaT = 2;           % 2K temperature difference

% calculation of Heat Transfer loop


for i= 1:numel(d_i)

switch true
    
    case Re(i) > 10000  
    Nu(i) = 0.024*Re(i)^0.8*Pr^(1/3);
    
    case 10000 > Re(i) > 2300
    Xi(i) = (1.8*log10(Re(i))-1.5)^-2; 
    Nu(i) = (Xi(i)/8*Re(i)*Pr)/(1+12.7*(Xi(i)/8)^(0.5)*(Pr^(2/3)-1))*(1+(d_i(i)/L)^(2/3));
    
    case Re(i) < 2300
    X(i) = L/(d_i(i)*Re(i)*Pr);
    Nu_O(i) = 3.657/(tanh(2.264*X(i)^(1/3)+1.7*X(i)^(2/3)))+0.0499/X(i)*tanh(X(i));
    Nu(i) = Nu_O(i) * 1 /(tanh(2.43*Pr^(1/6)*X(i)^6));
    
end

alpha(i) = Nu(i) * lambda/d_i(i)     % W/m^K
Q(i) = alpha(i) * A_ca(i) * DeltaT   % W

end

plot(d_i, Q)
xlabel('Diameter (m)')
ylabel('Heat Transfer in W')
title('Heat Transfer Depending on Diameters')
    
    
    
    
        
                       
                       









