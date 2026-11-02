# Family-state isolation correction

Machine: **easteregg2**.

The initial batch generator built s3k, w2lvar, and syn6 in one process. The w2lvar builder intentionally leaves its `VARIANT_G` overlay active for its physical runner. Syn6's extracted chemical matrices were correct, but its `global_G` metadata inherited w2lvar's heel/toe input gains of 0.5 rather than the accepted syn6 defaults of zero. Its inherited `renshaw` flag was also true, although the actual compiled graph already correctly contained no Renshaw cells.

The generator now resets pristine source parameters before every extraction and runs each family in a separate Python process. Syn6's physical-input metadata and flags are corrected. All neuron counts, chemical matrices, exact saved bilateral projections, and 64-step saved neural-current response parity remain unchanged. Source models, builders, learned weights, and manuscripts were not edited.

The first reflex experiment reused the affected metadata and therefore ran an altered-contact-input condition. Its failed physical-prefix audit prevented promotion. All numerical traces are retained at `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait2392_reflex_demo_v1`, together with `SYN6_LEDGER_BEFORE_G_ISOLATION.json` and `SOURCE_IMPLEMENTATION_V1.py`. It must not be called the frozen accepted standing condition.

The corrected reflex experiment is registered at `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait2392_reflex_demo_v2`. It independently restores pristine source defaults, requires them to equal the corrected metadata, explicitly requires heel/toe gains of zero, and gates further pulse cases on exact physical-prefix reproduction of the archived biological baseline. No gain fitting is performed.

A separate metadata-only unit label was corrected: source time step `0.002` is stored as `dt_s`, not `dt_ms`. This does not change calculations or figures. Native runner `q` traces are already in degrees, whereas `qfull` remains generalized SI/radian coordinates; the reflex postprocessor uses these native units without a second angular conversion.
