%% Single Izhikevich neuron (regular spiking)

clear; clc;

% --- Parameters ---
a = 0.02;      % excitatory neuron
b = 0.20;
c = -65;
d = 8;

dt = 1;        % ms, time step
T  = 1000;     % ms, total time
nt = T/dt;     

I  = 10;       % constant input current (adjust if no spikes)

% --- Initial conditions ---
V = -65;       % membrane potential (mV)
u = b*V;       % recovery variable

V_trace = zeros(nt,1);  % storage array, records V at every timestep

% --- Time loop ---
for t = 1:nt
    % Izhikevich equations 
    dV = 0.04*V^2 + 5*V + 140 - u + I;
    du = a*(b*V - u);

    V = V + dt*dV;
    u = u + dt*du;

    % Spike reset
    if V >= 30
        V = c;
        u = u + d;   
    end

    V_trace(t) = V;
end

% --- Plot membrane potential ---
time = (0:nt-1)*dt;

figure;
set(gcf,'Color','k');              % black figure background
plot(time, V_trace, 'Color',[0 1 1], 'LineWidth',1.5);   % cyan line

ax = gca;
ax.Color = [0 0 0];                % black axes background
ax.XColor = [1 1 1];               % white axis labels
ax.YColor = [1 1 1];

xlabel('Time (ms)');
ylabel('Membrane potential V (mV)');
title('Single Izhikevich neuron','Color',[1 1 1]);
grid on;
