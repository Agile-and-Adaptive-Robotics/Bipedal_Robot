# Goal-3 timestep convergence ladder (2026-09-23)

Same recorded 20-muscle activation control (spinal_run_s3k.npz,
2 ms zero-order hold), plant integrated at 2 / 1 / 0.5 ms.
Shared initial state; divergence = integration error only.

## ground (contacts on)

- 2 ms vs 0.5 ms: overall RMS **0.439 deg**, max 11.897 deg
- 1 ms vs 0.5 ms: overall RMS **0.251 deg**, max 6.797 deg
- wall time: 2 ms 1s, 1 ms 2s, 0.5 ms 4s

## air (no ground contact)

- 2 ms vs 0.5 ms: overall RMS **0.191 deg**, max 3.945 deg
- 1 ms vs 0.5 ms: overall RMS **0.122 deg**, max 2.525 deg
- wall time: 2 ms 1s, 1 ms 2s, 0.5 ms 4s
