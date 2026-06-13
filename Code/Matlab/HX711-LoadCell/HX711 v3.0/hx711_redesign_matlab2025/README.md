# HX711 MATLAB 2025 App Redesign

## Files

- `HX711.m` - programmatic MATLAB App Designer-compatible class. You can run it with `app = HX711;` or paste/import into App Designer and save as `HX711.mlapp`.
- `basic_HX711.m` - Arduino add-on MATLAB library class.
- `src/basic_HX711.h` - Arduino add-on C++ header.

## Intended paths

Replace the original app source around:

`GitHub\Bipedal_Robot\Code\Matlab\HX711 v3.0\HX711 v3.0\HX711.mlapp`

with this source after backing up the original `.mlapp`.

The default save folder in the app is:

`GitHub\Bipedal_Robot\Testing_Data\2026_06_Festo\Flx_20mm\`

## Save format

The Save button writes:

`FlxTest$_#.mat`

where `$` is the series spinner and `#` is a 2-digit run number spinner, `00` through `99`.

Each MAT-file contains:

- `Data`: `750 x 6` numeric array.
  - Column 1: `Time_s`
  - Column 2: `RawHX711_counts`
  - Column 3: `Force_N`
  - Column 4: `Pressure_kPa`
  - Column 5: `PressureVoltage_V`
  - Column 6: `Force_SelectedUnit`
- `Stats`: table with rows `Mean`, `Median`, `Mode`, `Min`, `Max`, `StdDev` and columns `Force`, `Pressure`.
- `Metadata`: struct containing pins, angles, units, calibration coefficients, and valve logic.
- `ColumnNames`: cell array naming `Data` columns.

## Key defaults

- Clock pin: `D2`
- HX711 data pin: `D3`
- Pressure sensor pin: `A0`
- Valve pins: `D11`, `D6`
- Force display unit: Newtons
- Sample target: `750`

## Valve logic

- Increase pressure: `D11 High`, `D6 High`
- Maintain pressure: `D11 Low`, `D6 High`
- Decrease pressure / 0 kPa maintain: `D11 Low`, `D6 Low`

## Calibration changes

Load cell:

- You can still use `Tare` and `Scale Factor`.
- Or enter known zero offset/tare and calibration slope in `Known LC Cal`, then click `Apply Known Load-Cell Cal`.

Pressure sensor:

- Click `Run 7-Point Pressure Cal`.
- The app prompts for actual gauge kPa at nominal regulator setpoints: `0, 200, 300, 400, 500, 600, 620 kPa`.
- The app reads A0 and fits `Pressure_kPa = a*Voltage_V + b` using `polyfit`.
- Or enter known `a` and `b`, then click `Apply Known Pressure Cal`.

## Notes before committing

1. The updated `src/basic_HX711.h` sets `libName = "basicHX711/basic_HX711"` to match `LibraryName` in `basic_HX711.m`.
2. If your existing Arduino add-on folder structure differs, preserve your folder structure and only update the contents of these two HX711 add-on files.
3. A `.mlapp` file is a packaged App Designer file. This deliverable is a source `.m` app class. Open it in MATLAB R2025, verify layout, then save/export as `.mlapp` if you want the packaged app format.
