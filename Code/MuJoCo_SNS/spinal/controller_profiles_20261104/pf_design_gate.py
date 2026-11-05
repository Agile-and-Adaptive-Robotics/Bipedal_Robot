"""Explicit hold on downstream fitting of the known failed PF revisions.

Replay/plotting and isolated PF-repair experiments remain available. Clear this
hold only when a revised implemented PF design has declared targets, measured
timing acceptance, and fresh-process evidence. Historical fits are unchanged.
"""
def require_pf_design_ready(family):
 reasons={
  's3k':'Four PF channels collapse toward two patterns; numerical PF onset/duration targets are not declared.',
  'syn6':'Primitive timing audit failed; S5 remains active throughout the cycle and the other channels collapse toward shared E/F profiles.'}
 if family in reasons:
  raise RuntimeError('PF design hold: '+reasons[family]+' Repair and verify PF timing before MN/readout or body optimization. Historical replay is still allowed.')
