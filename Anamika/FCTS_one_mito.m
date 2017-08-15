%FTCS implementation of one-mitochondria sink diffusion equation
%Dimensionless groups : 
%Delta = half of the extent of mitochondria
%Where actual extent of mitochondria = 2 micron
%1/K_g = dimensionless wrt time
%At K-g = 10^5 glucose molecules per second, 1/K_g = 1/10^5 seconds
%C_0 = dimensionless parameters
%discretization:
%dx = 0.1, dt = 10^-4
%All other parameters are scaled accordingly
%Variables with subscript a are actual values from literature in terms of
%microns and seconds. All other values are scaled values for simulation
clear;
dx = 0.1; 
delta_a = 1; %microns
delta = 1; %dimensionless group
k_g_a = 10^5; %glucose molecules per second
k_g = 1; %dimensionless group
dt = 10^-6; 
C_0_a = 1.5 * 10^6; %number of molecules
C_0 = 1;
dimless_t = 1/(k_g_a/C_0_a);
L_a = 500; %microns
L = L_a;
v_a = 2; %microns per second
v = v_a / (dimless_t /delta_a) ;
D_a = 10; %micron square per second
D = D_a /((delta_a)^2/dimless_t);
k_sg = 1; %constant for glucose dependent stopping rate k_s
k_w = 1; %Glucose independent walking rate
t = 0.01; 
steps = int64(t/dt);
G = zeros(steps,L/dx);
space_steps = int64(L/dx);
G(1,:) = C_0 ; %Initializing glucose field
m_pos = zeros(1,steps); %mitochondria position at every time step
m_pos(1) = delta;
m_run = 1; m_stop = 0; %assume mito is running initially
%calculating mitochondria position:
for n = 1:1:steps
    G_index = floor(m_pos(n)/dx);
    Gs = G(n,G_index);
    if(isnan(Gs))
       break;
    end
    if(m_run)
        p = rand;
        if(p>(k_sg*Gs/k_w + k_sg*Gs))
            m_stop = 1;
            m_run = 0;
            m_pos(n+1) = m_pos(n);
        else
            m_pos(n+1) = m_pos(n)+ (v * dt);
        end
    else
        p = rand;
        if(p>(k_w/(k_sg*Gs + k_w)))
            m_stop = 0;
            m_run = 1;
            m_pos(n+1) = m_pos(n)+ (v*dt);
        else
            m_pos(n+1) = m_pos(n);
        end
    end
    sink = 0;
    G_dash_dash(1) = 0;
    G_dash_dash(space_steps) = 0;
    G(:,1) = C_0;
    for j = 2:1:(space_steps-1)
        %evaluate the presence of sink for position j
        if (m_pos(n)-delta<(j*dx)<m_pos(n)+delta)
            sink = 1;
        end
        G_dash_dash(j) = (G(n,j+1) - 2*G(n,j) + G(n,j-1))/((dx)^2);
        G(n+1,j) = ((D * (G_dash_dash(j))) - (k_g * G(n,j) * sink) *dt) + G(n,j) ; 
    end
end

 
    







