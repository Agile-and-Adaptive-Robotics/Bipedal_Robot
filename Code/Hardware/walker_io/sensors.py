"""SKELETON: sensor readers for the walker. Every read_* returns the raw signal
in SI (kPa, m, rad, rad/s); sign conventions live in ONE place per device so the
afferent encoder stays clean. Hardware access points marked TODO(hw).
"""
from __future__ import annotations
import struct, time
from protocol import now


class PressureSensor:
    """BPA line pressure, 0.5-4.5 V ratiometric transducer -> kPa."""

    def __init__(self, name: str, adc_ch: int, kpa_per_v: float = 172.5, v_off: float = 0.5):
        self.name, self.adc_ch = name, adc_ch
        self.kpa_per_v, self.v_off = kpa_per_v, v_off

    def read_kpa(self) -> float:
        v = self._read_voltage()
        return max(0.0, (v - self.v_off) * self.kpa_per_v)

    def _read_voltage(self) -> float:
        # TODO(hw): ADC (ADS1256/ADS8688 on SPI) channel read
        raise NotImplementedError("ADC bring-up")


class Encoder:
    """Joint encoder -> theta, dtheta. Quadrature on GPIO/counter, or CAN."""

    def __init__(self, joint: str, counts_per_rad: float, sign: float = 1.0,
                 theta0: float = 0.0):
        self.joint, self.cpr, self.sign, self.theta = joint, counts_per_rad, sign, theta0
        self._last_counts, self._last_t, self.dtheta = None, None, 0.0

    def read(self) -> tuple[float, float]:
        c, t = self._read_counts(), now()
        if self._last_counts is not None:
            dth = self.sign * (c - self._last_counts) / self.cpr
            dt = max(1e-4, t - self._last_t)
            self.dtheta = 0.8 * self.dtheta + 0.2 * (dth / dt)     # light LPF
            self.theta += dth
        self._last_counts, self._last_t = c, t
        return self.theta + self.theta0, self.dtheta

    def _read_counts(self) -> int:
        # TODO(hw): quadrature counter / CAN read
        raise NotImplementedError("encoder bring-up")


class LiquidWire:
    """Muscle length sensor (Micklam). adc -> meters via per-sensor poly fit."""

    def __init__(self, muscle: str, adc_ch: int, coeffs: tuple[float, ...]):
        self.muscle, self.adc_ch, self.coeffs = muscle, adc_ch, coeffs
        self.L = None

    def read_m(self) -> float:
        a = self._read_adc()
        L = sum(c * a ** n for n, c in enumerate(self.coeffs))      # adc -> m
        self.L = L if self.L is None else self.L + 0.1 * (L - self.L)
        return self.L

    def _read_adc(self) -> float:
        # TODO(hw): same ADC chain as pressure, dedicated channels
        raise NotImplementedError("liquid wire bring-up")


class Imu:
    """Pelvis/trunk IMU: quaternion + gyro in body frame."""

    def __init__(self, site: str):
        self.site = site

    def read(self) -> dict:
        # TODO(hw): SPI/I2C IMU driver + fusion ( onboard or complementary filter)
        raise NotImplementedError("imu bring-up")


class InsoleArray:
    """Capacitive insole: heel/toe contact booleans first, CoP later.
    Front-end MCU streams over USB; this class parses the stream."""

    def __init__(self, side: str):
        self.side = side
        self.contact = {"heel": 0, "toe": 0}

    def read_contact(self) -> dict[str, int]:
        # TODO(hw): parse MCU stream; thresholds from the bench calibration
        raise NotImplementedError("insole bring-up")
