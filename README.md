# Honors Thesis: Bayesian Rank Aggregation and Distributed Computing

**A Statistical Study of Byzantine Consensus via Bayesian Rank Aggregation**

Departmental Honors in Applied Mathematics, Davidson College, 2025-26

- **Student:** Chadi Bsila
- **Advisors:** Hammurabi Mendes and Kaiwen (Kevin) Wang

## Overview

In distributed systems, a group of processes must agree on a shared outcome even when some of them are Byzantine (faulty or adversarial). Classical approaches to Byzantine consensus are combinatorial and rely on worst-case, fault-tolerant protocols that can be overly pessimistic. This project brings a statistical method to that setting: each process holds a noisy ranking of a set of items, and processes aggregate rankings using Bayesian rank aggregation as implemented in the [BiGER](https://github.com/kevin931/BiGER) R package.

This repository contains the simulation code used to study how well the aggregated rankings recover a ground-truth ranking as the share of Byzantine processes grows, how different communication structures affect robustness, and how quickly good processes converge over multiple rounds of communication.

The work has implications for voting theory and blockchains, where many parties must agree on a ranking or ordering in the presence of faulty or malicious participants.

## Repository contents

| File | Description |
|---|---|
| `algorithms.R` | Algorithm and network definitions: ranking simulation, single-round and multi-round aggregation, the intersection models (A, A', B, B', C, C'), and the Kendall's tau and Spearman evaluation metrics. |
| `simulations.R` | Main simulation script: one-round performance across good/Byzantine mixtures, convergence diagnostics (Gelman-Rubin), and sensitivity to the correlation parameter rho. Runs in parallel and writes CSV results to `simulation_outputs/`. |
| `multi_round_simulations.R` | Multi-round experiment script: runs repeated simulations over intersection models and values of `k`, aggregates mean Kendall's tau (or Spearman) per round with confidence intervals, plots the results, and computes the Gelman-Rubin diagnostic on process trajectories. |

## Reproducibility

Every simulation sets its own random seed, so results are reproducible for a given R and package version. Parallel runs use `future.seed = TRUE` for independent random-number streams.

## Results

All main experiments use n = 300 processes with t = 99 Byzantine processes (t < n/3) and m = 50 items. They are run under a worst-case, bimodal setting in which Byzantine processes collude around a single competing ranking. Each configuration is repeated 50 times.

- **One round of BiGER stays robust under heavy Byzantine presence.** As the number of Byzantine inputs seen by each good process grows from 1 to 99, the average Kendall's tau against the ground truth falls from 0.934 to 0.495, and the average Spearman correlation falls from 0.991 to 0.673. The output degrades gradually rather than collapsing, even when faulty processes make up nearly half of the inputs.
- **Posterior uncertainty identifies Byzantine processes.** BiGER's posterior uncertainty is clearly higher for Byzantine rankings than for good ones, and the resulting trust matrix separates the two groups after a single round.
- **Communication structure matters.** Intersection model A', which forces Byzantine processes into every good process's view, is by far the most adversarial. Tighter, gossip-like models (B and C) perform best, even compared with the random baseline A, and C' is much less adversarial than A'. Enforcing shared information between good processes limits the influence of Byzantine outputs.
- **Multiple rounds keep Byzantine influence in check.** In the multi-round procedure, each process accumulates the rankings it receives and uses BiGER's posterior uncertainty as trust scores, so rankings from Byzantine processes are down-weighted or ignored in later rounds. This was simulated under all six intersection models.
- **Good processes stabilize within a few rounds.** Using the Gelman-Rubin statistic (PSRF) as a convergence heuristic over 10 rounds and 50 repetitions, mean PSRF is between 0.991 and 1.010 for every intersection model, well within the usual 1.0 to 1.2 range for consensus.

## Acknowledgements

Thanks to my thesis advisors, Hammurabi Mendes and Kaiwen (Kevin) Wang, and to the Davidson College Department of Mathematics and Computer Science. This work builds on the [BiGER](https://github.com/kevin931/BiGER) package.
