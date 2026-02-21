clear; clc;

% Network setup 
N  = 100;            % total number of neurons
Ne = round(0.8*N);   % number of excitatory neurons
Ni = round(0.2*N);   % number of inhibitory neurons

dt = 1;              % simulation time step (ms)
T  = 1000;           % total simulation time (ms)
nt = T/dt;           % number of simulation time steps

% Izhikevich parameters for cortical E/I neurons
a = [0.02*ones(Ne,1); 0.10*ones(Ni,1)]; 
b = [0.20*ones(Ne,1); 0.20*ones(Ni,1)]; 
c = -65*ones(N,1);
d = [8*ones(Ne,1); 2*ones(Ni,1)];

V = -65*ones(N,1);   % initial membrane potential
u = b.*V;            % initial recovery variable 

% FIXED: Per-neuron synaptic weights (N x N matrix)
S = zeros(N,N);                     % Full connectivity matrix
S(:,1:Ne) = 0.5*rand(N,Ne);         % Excitatory → all positive (Dale's law)
S(:,Ne+1:N) = -rand(N,Ni);          % Inhibitory → all negative (Dale's law)

% Synapse parameters (Humphries 2009)
hE = zeros(N,1);        % excitatory gating variable (per postsynaptic neuron)
hI = zeros(N,1);        % inhibitory gating variable (per postsynaptic neuron)
tauE = 5;               % AMPA-like decay (ms)
tauI = 10;              % GABA-like decay (ms)
EE = 0;                 % excitatory reversal potential (mV)
EI = -70;               % inhibitory reversal potential (mV)
gE = 0.10;              % excitatory conductance 
gI = 0.10;              % inhibitory conductance

firings = []; % stores all spikes [time, neuron_id]

%% Time loop 
for t = 1:nt
    % External Poisson-like noise (higher for E neurons)
    I_ext = [5*randn(Ne,1); 2*randn(Ni,1)]; 
    
    % FIXED: TRUE PER-NEURON SYNAPTIC CURRENTS
    % Each postsynaptic neuron i gets: I_syn(i) = sum_j S(i,j) * spike_j(t)
    Sspike = zeros(N,1);                    % spike vector this timestep
    IE = gE * hE .* (EE - V);               % excitatory synaptic current (per neuron)
    II = gI * hI .* (EI - V);               % inhibitory synaptic current (per neuron)
    I_syn = IE + II;                        % total synaptic current (per neuron)
    I = I_ext + I_syn;                      % total current (per neuron)
    
    % Izhikevich model integration
    dV = 0.04*V.^2 + 5*V + 140 - u + I;
    du = a.*(b.*V - u);
    V = V + dt*dV;
    u = u + dt*du;

    % Spike detection
    fired = find(V >= 30);
    
    % Update synaptic gating variables (per postsynaptic neuron)
    if ~isempty(fired)
        Sspike(fired) = 1;                      % presynaptic spike indicators
        % PER-NEURON excitatory input: SE(i) = sum_{j in E} S(i,j)*Sspike(j)
        SE = S(:,1:Ne) * Sspike(1:Ne);          % N-vector of E inputs
        SI = S(:,Ne+1:end) * Sspike(Ne+1:end);  % N-vector of I inputs
        
        % Gating dynamics (Humphries 2009) - now truly per postsynaptic neuron
        hE = hE + dt * (-hE/tauE + SE);         % AMPA gating
        hI = hI + dt * (-hI/tauI + SI);         % GABA gating
        
        % Store spikes and reset
        firings = [firings; t*ones(length(fired),1), fired]; 
        V(fired) = c(fired);
        u(fired) = u(fired) + d(fired);
    else
        % Still update gating variables even without spikes (decay only)
        hE = hE + dt * (-hE/tauE);
        hI = hI + dt * (-hI/tauI);
    end
end

%% Raster plot
figure; set(gcf,'Color','k');
plot(firings(:,1), firings(:,2), '.', 'Color',[0 1 1], 'MarkerSize', 2);
ylim([0 N]); 
ax = gca; ax.Color = [0 0 0]; ax.XColor = [1 1 1]; ax.YColor = [1 1 1];
xlabel('Time (ms)'); ylabel('Neuron index');
title('Izhikevich E/I Network: Per-Neuron Synaptic Currents','Color',[1 1 1]);
