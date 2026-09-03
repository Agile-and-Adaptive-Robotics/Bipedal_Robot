# Deferred human / robot frame work

Recorded from the user's instructions. This update does NOT implement these
changes or alter transformation construction in any model/result script.

## Preserve explicit transformation calculations

Do not remove, relocate, or replace the user's explicit transformation-matrix
construction with a context-loading shortcut in other scripts. Ask before any
such refactoring. Human and robot ICR translations use different functions;
their knee rotation matrices are the same. Keep those frames distinct.

## Robot tibia points to human ICR frame

The user supplied this translation and direction:

```matlab
Htot1 = [-16.42, -47.11, 0]/1000; % metres
T_H_t1(:,:,i) = RpToTrans(eye(3), Htot1');
```

`T_H_t1(:,:,i)` is to translate robot tibia/t1 points into the human ICR frame.
Do not silently reverse its sign or substitute the robot's angle-dependent
t1-to-ICR transform. Confirm its exact place in the model's transform chain
when implementing that later task. It is not a camera or display transform.

## Current axes decision

The supplied knee rotation is Rz(phi): motion occurs in native X-Y and leaves
native Z invariant. Keep native X, Y and Z in plots. Both `DisplayRotation` and
`DisplayAxisMap` are identity for the extensor calls. `CameraOrbitDeg=90` changes
the view about Y without changing plotted coordinates. Do not reintroduce the
automatic `[-z,y,-x]` mapping as a supposed anatomical correction.
