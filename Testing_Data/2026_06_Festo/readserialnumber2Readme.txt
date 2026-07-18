With the current MATLAB code, pressing O and then V again starts a fresh unsaved test using the same next file number.

The behavior is:

V clears recordedData, turns the valves on, and begins a new recording.
O turns the valves off, but currently does not stop recording or clear the buffer.
S is the only command that selects a filename and saves it.

Therefore, an unsaved setup sequence such as:

V → observe → O → adjust setup → V → observe → O

does not consume any test numbers. When you eventually press S, MATLAB will still use the next available number.

One detail: after you press O, MATLAB continues storing valve-off data until you press V, S, or quit. Pressing V clears all of that unsaved data anyway.

Make O explicitly cancel an unsaved setup test

Because your intended workflow is S before O for a real test, it may be clearer to make O stop and discard any unsaved recording. Replace your current "o" case with:

case "o"
    writeline(s, "O");

    % Stop and discard the current unsaved setup test
    recording = false;
    recordingPending = false;
    recordedData = zeros(0, expectedNumCols);

    fprintf([ ...
        "Sent O: requesting both valve outputs LOW.\n" ...
        "Current unsaved recording was discarded.\n"]);

Do not clear data there. If you already pressed S, data contains the saved test and the files have already been written.