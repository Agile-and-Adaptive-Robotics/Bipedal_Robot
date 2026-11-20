# Connectome editor

Install Python 3.10 and the packages in `requirements.txt` in a virtual environment. From this directory run `python app.py`, or activate the environment and open `Start_Connectome_Editor.cmd` on Windows. Use **File → Open controller** to load Syn6, the s3k and w2lvar standing configurations, or the native W2L108 model. Each uses the input files included in this repository.

Double-click a module to see the neurons it contains. Select a connection to inspect its parameters. **Save As** creates a working copy of a template. Native W2L projects use an AnimatLab project and its paired simulation file in the same directory. The editor retains native cell and connection records when saving these models.

To run the Python checks, use `python -m unittest discover -s tests -p "test_*.py"`. The canvas checks use Node.js and run with `node tests/test_editor.cjs`.
