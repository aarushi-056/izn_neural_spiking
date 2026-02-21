# Izhikevich E/I Neural Network Simulations

MATLAB simulations of a 100-neuron spiking network using the Izhikevich (2003) neuron model with conductance-based AMPA/GABA synapses. Built as part of a university computational neuroscience assignment.

---

## Overview

This project models the spiking activity of a biologically realistic network of 100 neurons — 80 excitatory and 20 inhibitory — using the Izhikevich (2003) simple spiking neuron model. The project progresses from a single isolated neuron through to a fully connected network with conductance-based synaptic dynamics, exploring how synaptic parameters shape network-level synchrony.

---

## Files

| File | Description |
|---|---|
| `iznp1.m` | Baseline network — all neurons receive the same constant drive, no synapses |
| `izn_synaptic_weights.m` | Adds weight matrix `S`; synaptic input via direct spike-weight product |
| `izn_100_steps.m` | Full AMPA/GABA conductance model with gating variables `hE`, `hI` |
| `izn_per_neuron.m` | Corrected per-neuron synaptic currents with proper Dale's law connectivity |
| `izn_subplots.m` | 6-panel figure sweeping excitatory/inhibitory conductance (`gE`/`gI`) |
| `izn_tau_subplots.m` | 6-panel figure sweeping synaptic time constants (`tauE`/`tauI`) |

---

## Model Details

### Neuron Model — Izhikevich (2003)

```
dV/dt = 0.04V² + 5V + 140 - u + I
du/dt = a(bV - u)
```

If V ≥ 30 mV: V → c, u → u + d

The parameters (a, b, c, d) were chosen to approximate the firing behaviour of biological cortical neurons. The membrane potential V is initialised at −65 mV (resting potential), and the recovery variable u is initialised as b·V. The spike reset at 30 mV implements the refractory period observed in real neurons — V resets to c (threshold), and u jumps by d to model the after-hyperpolarisation that follows a spike.

| Neuron type | a | b | c | d |
|---|---|---|---|---|
| Excitatory | 0.02 | 0.2 | −65 | 8 |
| Inhibitory | 0.10 | 0.2 | −65 | 2 |

### Synaptic Weight Matrix

```matlab
S = [0.5*rand(N,Ne), -rand(N,Ni)]
```

Excitatory neurons are assigned positive weights (0 to 0.5), mimicking AMPA/NMDA receptor conductances. Inhibitory neurons are assigned negative weights (−1 to 0), mimicking GABA receptor conductances. This implements Dale's law — neurons are either purely excitatory or purely inhibitory. Each postsynaptic neuron receives input proportional to the weighted sum of presynaptic spikes, so every neuron in the network receives a different effective stimulus depending on its connectivity — making the model more biologically realistic than a homogeneous constant-drive network.

### Conductance-Based Synaptic Model — Humphries et al. (2009)

```
I_syn = gE·hE·(EE − V)  +  gI·hI·(EI − V)
dhE/dt = (−hE + SE) / tauE
dhI/dt = (−hI + SI) / tauI
```

| Parameter | Value | Biological meaning |
|---|---|---|
| `tauE` | 5 ms | AMPA decay time constant |
| `tauI` | 10 ms | GABA decay time constant |
| `EE` | 0 mV | Excitatory reversal (AMPA permeable to Na⁺/K⁺) |
| `EI` | −70 mV | Inhibitory reversal (GABA-A selective for Cl⁻) |
| `gE` | 0.10–0.15 | AMPA/glutamate conductance strength |
| `gI` | 0.10–0.40 | GABA conductance strength |

`hE` and `hI` represent the total fraction of open AMPA and GABA-A receptor channels on each postsynaptic neuron respectively. The external drive uses Gaussian noise (`5·randn` for excitatory, `2·randn` for inhibitory neurons), reflecting the larger dendritic arbors and higher background synaptic bombardment of excitatory neurons in cortex.

---

## Key Findings

### Effect of gE / gI Balance on Synchrony

Increasing `gE` amplifies recurrent excitatory feedback, producing tighter synchrony and larger population bursts. Increasing `gI` suppresses burst size and resets network activity rhythmically. At high values of both (gE=0.15, gI=0.40), strong inhibition disrupts the excitatory rhythm, producing messier, less coherent bursting. This demonstrates that **peak network synchrony requires a balanced E/I ratio** — consistent with the neuroscience literature on cortical oscillatory dynamics.

### Effect of Synaptic Time Constants on Synchrony

Shorter time constants (faster synapses) generally produce tighter synchrony, since postsynaptic neurons integrate and respond to presynaptic spikes more rapidly. Very short values (tauE=1, tauI=2) produce near-total, almost epileptic synchrony. However, the relationship is nonlinear — halving both time constants at a slower baseline does not reproduce the synchrony seen at faster absolute timescales, indicating that **absolute values of tauE/tauI matter, not just their ratio**.

The most biologically realistic activity pattern (tauE=5, tauI=10) maintains moderate synchrony alongside asynchronous irregular firing. This balance matters: highly synchronous networks cannot perform flexible population coding, while fully asynchronous networks fail to propagate any signal. Disruption of E/I balance and synchrony has been implicated in Parkinson's disease, Alzheimer's disease, epilepsy, and autism (Uhlhaas & Singer, 2006).

---

## Requirements

- MATLAB (no additional toolboxes required)

---

## Usage

Each `.m` file is self-contained — open any script and run it directly to produce a raster plot. To explore parameter effects, run `izn_subplots.m` (conductance sweep) or `izn_tau_subplots.m` (time constant sweep).

---

## References

- Izhikevich, E. M. (2003). Simple model of spiking neurons. *IEEE Transactions on Neural Networks*, 14(6), 1569–1572.
- Humphries, M. D. et al. (2009). Dopaminergic control of the exploration-exploitation trade-off via the Basal Ganglia. *Frontiers in Neuroscience*.
- Dayan, P., & Abbott, L. F. (2001). *Theoretical Neuroscience*. MIT Press.
- Uhlhaas, P. J., & Singer, W. (2006). Neural synchrony in brain disorders: Relevance for cognitive dysfunctions and pathophysiology. *Neuron*, 52(1), 155–168.
