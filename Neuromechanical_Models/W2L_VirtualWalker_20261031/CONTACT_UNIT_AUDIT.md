# Contact-current unit verification

Session 1 subagent; easteregg2; 31 October 2026. This audit distinguishes a **specific polynomial adapter** from the user's screenshot of a **sigmoid with Amplitude=10 n**. The latter is a legitimate nano prefix; it must not be described as 10 A merely because its displayed magnitude is 10.

## Exact source fields

The four contact polynomial adapters in the named NEWER editable project contain:

```xml
<C Value="1" Scale="None" Actual="1"/>
<D Value="0" Scale="None" Actual="0"/>
```

Their source `DataTypeID` is ContactCount and destination `TargetDataTypeID` is ExternalCurrent. The supplied NEWER standalone instead has:

```xml
<SourceDataType>ContactCount</SourceDataType>
<TargetDataType>ExternalCurrent</TargetDataType>
<RobotIOScale>1</RobotIOScale>
<DelayBufferMode>0</DelayBufferMode>
<Gain>
  <Type>Polynomial</Type><UseLimits>False</UseLimits>
  <A>0</A><B>0</B><C>20</C><D>-4e-010</D>
</Gain>
```

Project adapter GUIDs differ from the supplied standalone's cafe0211–cafe0214 GUIDs. They are paired by physical source and target neuron, not by assuming identical adapter IDs. The source GUI has C Actual=1, not a hidden Actual=1e−9. No assertion of original author intent is made from these values.

## Official source-code chain

Public source used for audit:

- [Polynomial GUI gain](https://github.com/NeuroRoboticTech/AnimatLabPublicSource/blob/master/Libraries/AnimatGUI/DataObjects/Gains/Polynomial.vb), SaveSimulationXml lines 230–238, delegates each coefficient to its ScaledNumber serializer. Live parameter updates likewise use coefficient ActualValue. GUI coefficient magnitudes and explicit prefixes therefore must be distinguished.
- [Sigmoid GUI gain](https://github.com/NeuroRoboticTech/AnimatLabPublicSource/blob/master/Libraries/AnimatGUI/DataObjects/Gains/Sigmoid.vb), lines 194–197 and 252–260, uses the same actual-valued serialization/update mechanism. Its displayed Amplitude is the B parameter. A nano-scaled amplitude is not the unscaled polynomial C found above.
- [Native polynomial](https://github.com/NeuroRoboticTech/AnimatLabPublicSource/blob/master/Libraries/AnimatSim/PolynomialGain.cpp), lines 149–153 and 210–213: reads the scalar coefficients and evaluates the polynomial. UseLimits is false in this source, so output limiting is not providing a nano conversion.
- [Native adapter](https://github.com/NeuroRoboticTech/AnimatLabPublicSource/blob/master/Libraries/AnimatSim/Adapter.cpp), lines 819–822 and 844–872: obtains gain output, applies RobotIOScale (1 here), and forwards it to the target. DelayBufferMode=0 adds no alternate scaling. GetDataPointer exposes CalculatedVal and NextVal for direct measurement.
- [Integrate-and-fire neuron](https://github.com/NeuroRoboticTech/AnimatLabPublicSource/blob/master/Libraries/IntegrateFireSim/Neuron.cpp), lines 1283–1290: `m_fltAdapterI += (fltInput*1e9);` and `m_fltAdapterMemoryI = (m_fltAdapterI * 1e-9);`. The AdapterCurrent getter exposes that second value, in A, while the internal neuron calculations use nA.

Public source is not automatically proof of the installed DLL version. The installed-runtime measurement below closes that gap for this specific model.

## Direct installed-native measurement

Created a **separate copy** of the original supplied NEWER standalone. Its coefficients, neurons, and graph were unchanged. Only the output directory, 0.51 s simulation end time, and diagnostic chart were changed. For each contact adapter the chart records physical ContactCount, adapter CalculatedVal, NextVal, target-neuron AdapterCurrent, and target-neuron membrane voltage.

| Signal condition | CalculatedVal | NextVal | Neuron AdapterCurrent |
|---|---:|---:|---:|
| Contact=0 | −4e−10 | −4e−10 | −4e−10 A |
| Left toe contact=1 at 0.036 s | 20 | 20 | **20 A** |
| Right toe contact=1 at 0.063 s | 20 | 20 | **20 A** |

Native exit is 0. Both heel sensors remain at zero contact in this short interval; the two toe channels directly demonstrate the actual on-contact value. Thus there is no hidden nano multiplier along the measured runtime path. A bounded membrane-voltage plot alone does not establish reasonable input-current magnitude.

Evidence is in `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\build\native_contact_unit_probe_20261031\PROBE_REPORT.json`, with full raw chart, probe model, and NPZ. The builder is tools/probe_native_contact_units.py. Original source files were not edited by the probe.

## Candidate interpretations, not recovered intent

- **GUI1nA:** C=1 nA/contact, D=0. This preserves GUI scalar magnitudes under an explicit nA-unit reinterpretation.
- **Existing SIContact engineering candidate:** C=20 nA/contact, D=−0.4 nA. This follows the older runtime scalar 20 and offset choice. Relative to GUI1nA it changes gain by 20× and adds a baseline offset. It is **not** a purely numerical-unit-preserving GUI conversion.

Neither gain was specified by Ben or identified as an optimum. Both candidates retain the source editable graph and its left/right RG-extensor tonic drive **2/0 nA**, so neither is a symmetric-drive controller. Tested 20 nA artifacts/results are retained unchanged. No parameter overwrite occurred after the user's scale challenge; the subsequent probe added telemetry to a separate unmodified-source copy.
