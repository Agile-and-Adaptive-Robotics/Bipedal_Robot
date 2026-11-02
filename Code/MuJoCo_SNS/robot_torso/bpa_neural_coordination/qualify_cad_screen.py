"""Clarify the historical geometric screen without changing its numeric results."""
from pathlib import Path
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031')
def main():
 path=ROOT/'cad_bpa/EXPERIMENT_RECORD.md';text=path.read_text(encoding='utf-8');heading='## Series-compliance qualification'
 if heading not in text:
  text+='\n'+heading+'\n\nThe reported 53/54 and 49/54 length-envelope counts assume zero series extension. They describe a geometric reference screen, not a hard impossibility bound for a compliant actuator. A total route may exceed sleeve rest plus fixed fittings while the cable stretches and the sleeve remains inside its characterized contraction range. Pressure capacity and series stiffness must be included before accepting or rejecting a complete compliant design. Numeric screen outputs are unchanged. The corresponding Gait pressure-interface regression is `'+str(ROOT/'gait_bpa/series_domain_regression.json')+'`. Mounted CAD parameters remain incomplete; no CAD pressure-standing result is claimed.\n'
  path.write_text(text,encoding='utf-8')
 print(path)
if __name__=='__main__':main()
