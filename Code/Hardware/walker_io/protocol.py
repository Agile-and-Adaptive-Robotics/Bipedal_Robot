"""SKELETON: wire contract between the walker's Orin Nano and the lab PC.
UDP, JSON lines. See Documentation/Program_Workflow/BPA_Walker_Program.md section 6.
"""
from __future__ import annotations
import json, time
from dataclasses import dataclass, field, asdict
from typing import Any

ORIN_IP, LAB_IP = "192.168.10.1", "192.168.10.2"
TELEM_PORT, CMD_PORT = 50100, 50101
TELEM_HZ = 50
HEARTBEAT_HZ = 10
SENSOR_STALE_S = 0.020          # >20 ms without a sensor stream -> vent-all


@dataclass
class Telemetry:
    t: float
    phase: float = 0.0
    activations: dict[str, float] = field(default_factory=dict)   # pool -> a in [0,1]
    afferents: dict[str, float] = field(default_factory=dict)     # "Ia_VASTI_R" -> rate
    pressures: dict[str, float] = field(default_factory=dict)     # BPA name -> kPa
    joints: dict[str, list] = field(default_factory=dict)         # joint -> [th, dth]
    imu: dict[str, list] = field(default_factory=dict)            # imu -> [qx,qy,qz,qw,wx,wy,wz]
    contact: dict[str, int] = field(default_factory=dict)         # "heel_R" -> 0/1
    boom_load: float = 0.0

    def to_json(self) -> bytes:
        return json.dumps({"k": "telem", "d": asdict(self)}).encode()


@dataclass
class Command:
    cmd: str
    payload: dict[str, Any] = field(default_factory=dict)

    @staticmethod
    def override_gains(g: dict[str, float]) -> "Command":      # neural deletion
        return Command("override", {"G": g})

    @staticmethod
    def stim(neuron: str, I_nA: float, dur_s: float, t_s: float | None = None) -> "Command":
        return Command("stim", {"neuron": neuron, "I": I_nA, "dur": dur_s, "t": t_s})

    @staticmethod
    def valve_override(name: str, duty: float) -> "Command":   # manual valve
        return Command("valve_override", {"name": name, "duty": duty})

    @staticmethod
    def stop() -> "Command":
        return Command("stop", {})

    def to_json(self) -> bytes:
        return json.dumps({"k": "cmd", "c": self.cmd, "d": self.payload}).encode()


def decode(raw: bytes) -> dict:
    return json.loads(raw.decode())


class UdpLink:
    """Tiny UDP wrapper; lab PC and Orin each run one in each direction."""

    def __init__(self, local_port: int, remote: tuple[str, int] | None = None):
        import socket
        self.s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.s.bind(("0.0.0.0", local_port))
        self.remote = remote

    def send(self, payload: bytes) -> None:
        if self.remote:
            self.s.sendto(payload, self.remote)

    def recv(self, timeout: float = 0.0) -> dict | None:
        self.s.settimeout(timeout)
        try:
            raw, _ = self.s.recvfrom(65535)
            return decode(raw)
        except TimeoutError:
            return None


def now() -> float:
    return time.monotonic()
