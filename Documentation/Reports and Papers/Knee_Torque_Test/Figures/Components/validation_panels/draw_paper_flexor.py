"""Rebuild Figure 12(B) using the same formatting as Figure 12(A)."""

import argparse
from pathlib import Path
import sys

from scipy.io import loadmat

PLOT_SOURCE = Path(__file__).resolve().parent.parent / 'torque_plots'
sys.path.insert(0, str(PLOT_SOURCE))
from draw_plots import configure, validation


def draw(output, fonts):
    configure(fonts)
    data = loadmat(PLOT_SOURCE / 'plot_data.mat', variable_names=['props'], simplify_cells=True)
    validation(output, 'flexor', data['props'])


if __name__ == '__main__':
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument('output', type=Path)
    args.add_argument('--fonts', type=Path, default=Path('C:/Windows/Fonts'))
    opts = args.parse_args()
    draw(opts.output, opts.fonts)
