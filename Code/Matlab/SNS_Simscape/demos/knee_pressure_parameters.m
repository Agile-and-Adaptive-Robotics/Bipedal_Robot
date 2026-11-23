function p = knee_pressure_parameters()
% The 20 mm comparison specimen is the 450 mm group of Chapter 3 and
% Testing_Data/CallModelCompare.m. Both fixed-radius routes use this same
% measured sleeve: unpressurized clamp spacing and zero-load 620 kPa spacing.
p.restingLength = 0.450;
p.freeLength620 = 0.335;
p.referenceLength = (p.restingLength + p.freeLength620)/2;
p.referenceAngle = pi/4;
p.pressureSupply = 620;             % kPa, gauge
p.r_arm = 0.035;                    % m, equal and opposite constant arms
p.I_knee = 0.06;                    % kg m^2
p.b_knee = 0.5;                     % N m s/rad
p.K_knee = 0.5;                     % N m/rad
p.th0 = 0.26;                       % rad, return spring's reference angle
p.ROM = pi/2;
p.Tload = 0.5;                      % N m, positive flexion
% Gait2392 Thelen2003 constants. These shape the control signal only.
p.activationTime = 0.010;
p.deactivationTime = 0.040;
p.motorThreshold = -45;             % mV
p.motorWidth = 10;                  % mV, full range of the bounded voltage ramp
p.connectionScale = 0.05;           % common sensory gain selected with the knee unloaded
p.feedbackGain = 1;                 % 0 disconnects all six sensory-to-motor pathways
p.tapAmplitude = 2;                 % nA, delivered to the extensor Ia neuron
p.tapOnset = 30;                    % s
p.tapDuration = 0.030;              % s
p.initialActivationFlex = 0.55;
% Choose descending currents to balance the prescribed operating point using
% the measured force surface. During simulation the two library blocks
% independently recalculate force at the actual pressure and length.
rho = (p.restingLength-p.referenceLength)/(p.restingLength-p.freeLength620);
f620 = 620*1.4877*atan(0.0248*620*(p.restingLength-0.0075));
forcePerActivation = f620*exp(-1.32087718059*rho^2);
passiveTorque = p.Tload-p.K_knee*(p.referenceAngle-p.th0);
p.initialActivationExt = p.initialActivationFlex + passiveTorque/(p.r_arm*forcePerActivation);
offset = f620*0.257852586017*(exp(-6.4766142989*rho)-1);
force = offset + forcePerActivation*[p.initialActivationExt,p.initialActivationFlex];
% Steady nonspiking Ia and Ib afferents and all three inputs to each MN.
% The spiking experiment uses the SAME mechanics, excitations, and weights
% it inherits from the baseline outside its replaced sensory/synapse layer.
sensoryRest = -52; gSensory = 0.4; gMN = 0.5;
iaCurrent = 6*[p.referenceAngle/p.ROM,1-p.referenceAngle/p.ROM];
ibCurrent = min(10,max(0,force.*[0.0200,0.0222]));
ia = min(1,max(0,(sensoryRest + iaCurrent/gSensory + 45)/10));
ib = min(1,max(0,(sensoryRest + ibCurrent/gSensory + 45)/10));
voltage = p.motorThreshold+p.motorWidth*[p.initialActivationExt,p.initialActivationFlex];
drive = gMN*(voltage - sensoryRest) + p.connectionScale*((0.18*ia).*voltage ...
      + (0.20*ib + 0.15*fliplr(ia)).*(voltage+72));
p.desc_ext = drive(1);
p.desc_flex = drive(2);
end
