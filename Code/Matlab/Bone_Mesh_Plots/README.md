# Native axes and optional full-skeleton animation

## Install

Replace `AnimateKneeBoneMuscle.m` and `TestKneeStaticSkeleton.m` in your
`Bone_Mesh_Plots` folder with the two files in this package, using their exact
names without numbered suffixes. Keep your current `MuscleBonePlotting.m` and
`TestKneeBonePlotting.m` from the preceding package.

Copy the contents of `Knee_Extensor_Plot_Calls.m` over only the two plotting
sections you pasted from `Knee_Extensor_20mm.m`. It uses the same variables,
frame order, colors and GIF settings. It does NOT rebuild transformations or
replace your complete result script.

After replacing the plotter, run `clear AnimateKneeBoneMuscle` and `rehash`,
then recreate the figures using the corrected plotting calls. No new
optimization is needed. Existing open figures do not update themselves.

## Fix both display matrices, not the camera

Your supplied rotation `[0 0 -1;0 1 0;-1 0 0]` maps a row vector `[x y z]` to
`[-z y -x]`: native X is placed on displayed Z, and native Z on displayed X.
The supplied native knee rotation instead acts in X-Y and leaves Z unchanged.
The anatomical plot therefore needs the native coordinates:

```matlab
'DisplayRotation',eye(3), ...
'DisplayAxisMap',eye(3), ...
'CameraOrbitDeg',90, ...
```

Both display matrices default to identity in this version. There is no hidden
automatic X/Z swap, including when `DisplayAxisMap=[]`. Explicit custom maps
are still accepted and applied exactly once. Keeping the old nonidentity
`DisplayRotation` in a caller will still swap coordinates, so replace that
line too. The camera orbit changes only the viewing direction.

The previous tests passed because they checked the supplied `[-z y -x]` map;
that did not establish its anatomical correctness. The revised tests verify
native XYZ preservation, including with a 90-degree camera orbit.

## Skeleton selection

Add this argument to the animated call (not the shared list, unless you want
it to control both figures):

```matlab
'FullSkeleton',true, ... % false for femur and tibia only
```

True shows femur, tibia, pelvis/hip, sacrum, lumbar spine, and foot. The foot
uses the human tibia transform at each animation frame. Pelvis/sacrum/spine
remain fixed. Extra meshes load once; fixed limits are computed over the
entire selected sweep, including the foot. False loads only femur/tibia.

When `FullSkeleton` is omitted, static remains full skeleton and animation
remains femur/tibia-only. The old `StaticFullSkeleton` option remains supported
for static calls; an explicitly supplied `FullSkeleton` overrides it.

Static padding remains 0.15 m per side. Existing animation +X/Z padding,
colors, human/robot route transforms, moving contacts and endpoint handling
are preserved. No changes to optimization or model transformation builders.

## Deferred frame work

See `FRAME_NOTES.md` for the user's supplied `Htot1` translation and instruction
to preserve explicit transformation construction in other scripts. These notes
are saved for later, not implemented by this plotting update.

## Validation

The three MATLAB files pass MISS_HIT syntax parsing. Independent numerical
checks cover XYZ preservation, the native X-Y flexion plane, 100-pose foot
transforms, distance preservation, and whole-sweep bounds. MATLAB graphics
tests have not been run here because MATLAB is unavailable.

To run the revised nine-test suite on your machine:

```matlab
clear AnimateKneeBoneMuscle TestKneeStaticSkeleton
rehash
results = runtests('TestKneeStaticSkeleton.m');
disp(table(results))
```

Besides the existing placement and route checks, the suite covers camera-only
rotation, full-skeleton animation, and equality of all static/animated bone
coordinates at the same zero-angle pose. Tests do not regenerate the actual
project figures; recreate those with the plotting section afterward.
