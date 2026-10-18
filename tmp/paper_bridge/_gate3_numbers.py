# Gate 3 supervisor check: verify 10 load-bearing manuscript numbers against
# the bridge results ledger + chapter digest + dissertation chapter of record.
import re, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
R = r"D:\GitHub\Bipedal_Robot"
MAIN = R + r"\Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex"
LED  = R + r"\tmp\paper_bridge\mining\results_ledger.md"
DIG  = R + r"\tmp\paper_bridge\mining\chapter_digest.md"
CH   = R + r"\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex"

def lines_with(path, pats, label):
    txt = open(path, encoding="utf-8", errors="replace").read().splitlines()
    print(f"--- {label} ---")
    for p in pats:
        hits = [i+1 for i, l in enumerate(txt) if re.search(p, l)]
        print(f"  /{p}/ -> lines {hits[:6]}")

# (claim, Main pats, Ledger pats, Digest pats, Chapter pats)
CHECKS = [
 ("1. flexor training-mean 8.12 -> 1.56 N*m",
  [r"8\.12", r"1\.56", r"2\.29"], [r"8\.122", r"1\.557"], [r"8\.122 . 1\.557"], [r"8\.12"]),
 ("2. adopted flexor chi +3.94 / 3.998e4 / 1.473e4",
  [r"\+3\.94", r"3\.998e4", r"1\.473e4"], [r"\+3\.94", r"3\.998e4", r"1\.473e4", r"3\.937370e-03"],
  [r"\+3\.94 mm, .1 = 3\.998e4"], [r"\+3\.94"]),
 ("3. adopted extensor chi -6.38 / 0.158 + bounds",
  [r"-6\.38", r"\b0\.158\b", r"qtyrange\{-20\}\{5\}"], [r"-6\.38", r"0\.158", r"-20, \+5"],
  [r"-6\.38 mm, .3 = 0\.158"], [r"-6\.38"]),
 ("4. initial set 11.9mm/102/10.6 + GoF 1.616/0.0958/3.2527 + ext -10.9/0.123",
  [r"11\.9", r"102", r"10\.6", r"1\.616", r"0\.0958", r"3\.2527", r"-10\.9", r"\b0\.123\b"],
  [r"0\.0119", r"1\.02E\+05", r"1\.06E\+04", r"1\.616", r"3\.2527", r"-0\.0109", r"\b0\.123\b"],
  [r"11\.9 mm", r"1\.616"], [r"0\.0119", r"1\.616", r"0\.123"]),
 ("5. direct ext 35.45-109.4 vs 20.25-50.32; +5.02% at -20.2; endpoints/rest/tendon",
  [r"35\.45", r"109\.4", r"20\.25", r"50\.32", r"\+5\.02", r"-20\.2", r"827\.7", r"36\.8"],
  [r"35\.4525", r"109\.3595", r"20\.248", r"50\.321", r"0\.0502493", r"827\.7"],
  [r"35\.45.109\.4", r"\+5\.02% at .20\.2"], [r"35\.45", r"827\.7"]),
 ("6. cam +30.83/+39.77 winner +31.30/+42.34; 269-507 kPa; 96.9mm; 43.1mm; 59.4-102.8 mean 72.4",
  [r"30\.83", r"39\.77", r"31\.30", r"42\.34", r"269", r"507", r"96\.9", r"43\.1", r"59\.4", r"102\.8", r"72\.4"],
  [r"30\.83", r"39\.77", r"31\.30", r"42\.34", r"269.507", r"96\.9", r"43\.1", r"59\.4", r"72\.4"],
  [r"\+30\.83%.+39\.77%", r"269.507 kPa", r"96\.9 mm"], [r"30\.83", r"96\.9", r"72\.4"]),
 ("7. flexor redesign -81.23..-3.864 vs -23.16..-2.461; margin 0.0501/5.01% at -107",
  [r"-81\.23", r"3\.864", r"-23\.16", r"2\.461", r"0\.0501", r"\+5\.01", r"-107"],
  [r"-81\.2343", r"-3\.8643", r"-23\.1587", r"-2\.4612", r"0\.0501"],
  [r"81\.23 to .3\.864", r"\+5\.01% at .107"], [r"81\.23", r"0\.0501"]),
 ("8. mono-pulley 85.74/85.56; 4.51 (5.6%); off-axis 0.803; objectives 74.47/0.08465; factors 3.51/3.70",
  [r"85\.74", r"85\.56", r"4\.51", r"5\.6", r"0\.803", r"74\.47", r"0\.08465", r"3\.51", r"3\.70"],
  [r"85\.7402", r"85\.5621", r"4\.506", r"0\.8034", r"74\.4733", r"0\.0846509"],
  [r"85\.74 N.+85\.56 N", r"4\.51 N", r"3\.51/3\.70"], [r"85\.74", r"74\.47"]),
 ("9. hardware points 1.163/7.849/23.437; loads 3.7/25.2/76; 388.5/388/592 kPa",
  [r"1\.163", r"7\.849", r"23\.437", r"3\.7;25\.2;76", r"388\.5;388;592"],
  [r"1\.1630 7\.8494 23\.4367", r"3\.7 / 25\.2 / 76 N", r"388\.5 / 388 / 592"],
  [r"1\.163 N", r"23\.437 N"], [r"1\.163", r"388\.5"]),
 ("10. bi-pulley 0.9968/0.9288/3.423/0.2423/101; clearances 19.142/15.558/7.024/4.260; 40.2 at 804 N; 25 mm",
  [r"0\.9968", r"0\.9288", r"3\.423", r"0\.2423", r"\b101\b", r"19\.142", r"15\.558", r"7\.024", r"4\.260", r"40\.2", r"\b804\b", r"\b25\b"],
  [r"0\.996794", r"0\.928794", r"3\.42299", r"19\.142", r"15\.558", r"7\.024", r"4\.260", r"40\.2", r"\b804\b"],
  [r"soleus 0\.9968", r"4\.260.19\.14", r"40\.2 N"], [r"0\.9968", r"4\.260"]),
]

for claim, mp, lp, dp, cp in CHECKS:
    print(f"\n===== {claim} =====")
    lines_with(MAIN, mp, "Main.tex")
    lines_with(LED, lp, "results_ledger.md")
    lines_with(DIG, dp, "chapter_digest.md")
    lines_with(CH, cp, "30-torque.tex")
