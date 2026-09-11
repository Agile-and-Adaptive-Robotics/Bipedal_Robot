"""Live neural/muscle strip-chart window that runs beside the MuJoCo viewer.

`runner.py --view` opens this scope next to the 3D window: top strip shows
spinal-cord potentials (DRIVE, POSTURE, RG-E/F per side, PF cells), bottom
strip shows the muscle activations they drive. Real-time pacing makes the
traces scroll at wall-clock speed.

Matplotlib interactive mode (plt.pause pumps the GUI event loop), so this
coexists with the MuJoCo viewer window without threads.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


class NeuroScope:
    """Preallocated ring buffers + Line2D redraws of DECIMATED data
    (plotting every 11k sample makes matplotlib the bottleneck; every
    ~5th sample is plenty at 500 Hz sim rate)."""

    REDRAW_EVERY = 0.5   # s of sim time between redraws

    def __init__(self, nsteps: int, dt: float, top_labels: list[str],
                 bot_labels: list[str]):
        self.dt = dt
        self.n = nsteps
        self.t = np.zeros(nsteps)
        self.top = np.zeros((nsteps, len(top_labels)))
        self.bot = np.zeros((nsteps, len(bot_labels)))
        self.k = 0

        self.fig, (ax_t, ax_b) = plt.subplots(
            2, 1, figsize=(9, 6), sharex=True, num="spinal cord - live")
        self.fig.subplots_adjust(left=0.09, right=0.98, top=0.93,
                                 bottom=0.08, hspace=0.12)
        self.lines_top = [ax_t.plot([], [], lw=1.0, label=l)[0]
                          for l in top_labels]
        self.lines_bot = [ax_b.plot([], [], lw=1.0, label=l)[0]
                          for l in bot_labels]
        ax_t.set_ylabel("neural [mV]")
        ax_b.set_ylabel("activation [0-1]")
        ax_b.set_xlabel("t [s]")
        for ax in (ax_t, ax_b):
            ax.grid(True, alpha=0.25)
            ax.set_xlim(0, max(1.0, nsteps * dt))
        ax_t.set_ylim(-4.5, 5.5)     # fixed: redraws stay cheap
        ax_b.set_ylim(-0.05, 1.05)
        ax_t.legend(fontsize=7, ncol=4, loc="upper left")
        ax_b.legend(fontsize=7, ncol=4, loc="upper left")
        ax_t.set_title("spinal network (live)", fontsize=10)
        plt.ion()
        plt.show(block=False)
        self._last_draw = -10.0

    def update(self, t: float, top_vals, bot_vals, force: bool = False):
        """Append one sample; redraw decimated data infrequently."""
        if self.k >= self.n:
            return
        self.t[self.k] = t
        self.top[self.k] = top_vals
        self.bot[self.k] = bot_vals
        self.k += 1
        if not force and t - self._last_draw < self.REDRAW_EVERY:
            return
        self._last_draw = t
        n = self.k
        s = max(1, n // 1200)     # decimate to <= ~1200 points per line
        for i, ln in enumerate(self.lines_top):
            ln.set_data(self.t[:n:s], self.top[:n:s, i])
        for i, ln in enumerate(self.lines_bot):
            ln.set_data(self.t[:n:s], self.bot[:n:s, i])
        try:
            self.fig.canvas.draw_idle()
            plt.pause(0.001)
        except Exception:
            pass  # window closed by the user; keep simulating

    def alive(self) -> bool:
        return plt.fignum_exists("spinal cord - live")

    def close(self):
        try:
            plt.ioff()
            plt.show(block=False)  # leave the finished trace up for review
        except Exception:
            pass
