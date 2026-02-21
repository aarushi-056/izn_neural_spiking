clear; clc; 

%Network setup 
N = 100; 
Ne = round(0.8*N);   %number of excitatory neurons
Ni = round(0.2*N);    %number of inhibitory neurons

dt = 1;              % simulation time step
T  = 1000;           % total simulation time
nt = T/dt;           %number of simulation time steps

% --- Izhikevich parameters from assignment ---
% Excitatory: a=0.02, b=0.2, c=-65, d=8
% Inhibitory: a=0.1,  b=0.2, c=-65, d=2
a = [0.02*ones(Ne,1); 0.10*ones(Ni,1)];
b = 0.20*ones(N,1);
c = -65*ones(N,1);
d = [8*ones(Ne,1); 2*ones(Ni,1)];

% --- Initial conditions ---
V = -65*ones(N,1);        % membrane potentials
u = b.*V;                 % recovery variables

% Connectivity: excitatory columns positive, inhibitory columns negative
S = [0.5*rand(N,Ne), -1*rand(N,Ni)];   % base weights 

% --- External current + synaptic input ---
I_ext = [5*randn(Ne,1); 2*randn(Ni,1)];  % noise 
I_syn = S * (V >= 30);                   % synaptic input from spikes 
I = I_ext + I_syn;                       % total current 

% --- Storage for spikes ---
firings = [];              % [time, neuron_index]

% --- Main simulation loop ---
for t = 1:nt
    
    % Update synaptic input each step (ADDED)
    I_syn = S * (V >= 30);
    I_ext = [5*randn(Ne,1); 2*randn(Ni,1)];  % fresh noise
    I = I_ext + I_syn;
    
    % Izhikevich equations
    dV = 0.04*V.^2 + 5*V + 140 - u + I;
    du = a.*(b.*V - u);
    V  = V + dt*dV;
    u  = u + dt*du;
    
    % Spike threshold and reset
    fired = find(V >= 30);
    if ~isempty(fired)
        firings = [firings; t*ones(length(fired),1), fired];
        V(fired) = c(fired);
        u(fired) = u(fired) + d(fired);
    end
end

% --- Raster plot ---
figure;
set(gcf,'Color','k');
if ~isempty(firings)
    plot(firings(:,1), firings(:,2), '.', 'Color',[0 1 1]);
    ylim([0 N]);  
end

xlabel('Time (ms)');
ylabel('Neuron index');
title('Izhikevich E-I network with synaptic weights');
