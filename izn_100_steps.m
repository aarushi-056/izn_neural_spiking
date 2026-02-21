%%Izhikevich network with simple AMPA/GABA synapses

clear; clc;

% Network setup 
N  = 100;            %total number of neurons
Ne = round(0.8*N);   %number of excitatory neurons
Ni = round(0.2*N);    %number of inhibitory neurons

dt = 1;              % simulation time step
T  = 1000;           % total simulation time
nt = T/dt;           %number of simulation time steps

% Izhikevich parameters
%cortical excitatory neuron parameters defined by Izhikevich
a = [0.02*ones(Ne,1); 0.10*ones(Ni,1)]; %vectors containing izhikevich 
b = [0.20*ones(Ne,1); 0.20*ones(Ni,1)]; %parameters for all neurons
c = -65*ones(N,1);
d = [8*ones(Ne,1); 2*ones(Ni,1)];

V = -65*ones(N,1); %membrane potential for each neuron
u = b.*V;          %recovery variable 

% Connectivity: excitatory columns positive, inhibitory columns negative
S = [0.5*rand(N,Ne), -1*rand(N,Ni)];   % base weights

% Synapse parameters (From Humphries paper)
hE = zeros(N,1);        % excitatory gating variable
hI = zeros(N,1);        % inhibitory gating variable

tauE = 5;               % ms, AMPA-like decay
tauI = 10;              % ms, GABA-like decay

EE = 0;                 % mV, excitatory reversal
EI = -70;               % mV, inhibitory reversal

gE = 0.15;              % excitatory conductance 
gI = 0.10;              % inhibitory conductance

% Spike storage
firings = []; %stores all spikes as rows 

%% Time loop 
for t = 1:nt %produces time evolution of the network

    % External noise current(random gaussian noise) %necessary step? why? 
    I_ext = [5*randn(Ne,1); 2*randn(Ni,1)]; 
    %excitatory neurons get slightly higher noise

    % Synaptic currents from previous hE, hI
    %gE maximal excitatory conductance, hE is gating variable
    II = gI .* hI .* (EI - V); %gI maximal inhibitory conductance, hI is gating variable
    IE = gE .* hE .* (EE - V);
    I  = I_ext + IE + II;      
    %total input current = external noise + excitatory synapses + inhibitory synapses

    % Izhikevich equations
    dV = 0.04*V.^2 + 5*V + 140 - u + I; %membrane potential
    du = a.*(b.*V - u);                 %recovery variable
    V  = V + dt*dV;
    u  = u + dt*du;

    % Detect spikes 
    fired = find(V >= 30); %spikes happen after crossing threshold 30 mV

    % Update synaptic gating based on these spikes
    Sspike = zeros(N,1); 
    Sspike(fired) = 1; %if neuron s(j) fires sspike = 1, or sspike =0

    SE = S(:,1:Ne)     * Sspike(1:Ne);      % excitatory input
    SI = S(:,Ne+1:end) * Sspike(Ne+1:end);  % inhibitory input
   %represent the terms in the equation 
   % (Theoretical Neuroscience 'Dayan and Abbot')
    hE = hE + dt*(-hE + SE)/tauE;
    hI = hI + dt*(-hI + SI)/tauI;
    % Store spikes and reset
    if ~isempty(fired)
        firings = [firings; t*ones(length(fired),1), fired]; 
        V(fired) = c(fired); % V is reset to c -> threshold
        u(fired) = u(fired) + d(fired);  %u is reset to u + d -> recovery 
    end
end

%% Raster plot 
figure;
set(gcf,'Color','k');
if ~isempty(firings)
    plot(firings(:,1), firings(:,2), '.', 'Color',[0 1 1]); % cyan 
end
ylim([0 N]); 
ax = gca;
ax.Color  = [0 0 0]; %because I'm using this in dark mode,setting axis colours 
ax.XColor = [1 1 1];
ax.YColor = [1 1 1];

xlabel('Time (ms)');
ylabel('Neuron index');
title('Izhikevich network with AMPA/GABA synapses','Color',[1 1 1]);
