"""Wall-clock limits for Ben's explicitly authorized nine-hour campaign."""
from datetime import datetime, timezone, timedelta
import math

START = datetime(2026, 11, 2, 10, 13, 22, tzinfo=timezone.utc)
DEADLINE = START + timedelta(hours=9)
TRAINING_STOP = DEADLINE - timedelta(minutes=90)
VERIFICATION_STOP = DEADLINE - timedelta(minutes=45)


def remaining(until=DEADLINE, now=None):
    return max(0., (until - (now or datetime.now(timezone.utc))).total_seconds())


def training_plan(pilot_wall_s, requested_trials=24, workers=1, now=None):
    """Use measured wall time, a 2x allowance, and a fixed small-run cap."""
    if not math.isfinite(pilot_wall_s) or pilot_wall_s <= 0:
        raise ValueError("A positive measured pilot wall time is required")
    if pilot_wall_s >= 900.:
        raise ValueError("Pilot exceeds the bounded per-trial limit; reduce the study before training")
    if workers not in (1, 2, 3) or not 1 <= requested_trials <= 24:
        raise ValueError("Small run: 1–24 trials and at most three workers")
    available = min(5400., max(0., remaining(TRAINING_STOP, now) - 120.))
    safe_trial_s = max(5., 2. * pilot_wall_s)
    trials = min(requested_trials, int(available * workers // safe_trial_s))
    return dict(max_trials=trials, workers=workers, wall_budget_s=available,
                measured_pilot_wall_s=pilot_wall_s, safe_trial_wall_s=safe_trial_s,
                per_trial_timeout_s=min(900., max(30., 3. * pilot_wall_s)),
                training_stop_utc=TRAINING_STOP.isoformat(),
                final_deadline_utc=DEADLINE.isoformat())


def require_time(estimated_wall_s, *, training=False, now=None):
    if not math.isfinite(estimated_wall_s) or estimated_wall_s <= 0:
        raise ValueError("A finite positive wall-time estimate is required")
    stop = TRAINING_STOP if training else VERIFICATION_STOP
    if estimated_wall_s > remaining(stop, now):
        raise TimeoutError("Job would consume the reserved verification/report time")
