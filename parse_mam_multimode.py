#!/usr/bin/env python3
"""
Parse the clean MAM no-deposition debug log and the STD HETP log,
produce a multi-mode NH4/NO3/SO4 budget table.
"""

import re, sys, math

MAM_LOG = "/cluster/work/climate/fsolmon/run_mosaic_v14_7/print_mam_nodep.log"
STD_LOG = "/cluster/work/climate/fsolmon/run_std_14.7/print_hetp_nodep.log"

# ──────────────────────────────────────────────────────────────────────────────
# Parser helpers
# ──────────────────────────────────────────────────────────────────────────────

def grab(line, key):
    """Return float value from a debug line containing 'key =' if key matches."""
    m = re.search(r'=\s*([0-9Ee.+\-]+)', line)
    return float(m.group(1)) if m else None

# Labels we care about in the MAM log
MAM_BEFORE_KEYS = {
    "RH":        "RH",
    "HNO3":      "HNO3_b",
    "NH3":        "NH3_b",
    "MAMSO41":   "SO41_b",
    "MAMNO31":   "NO31_b",
    "MAMNH41":   "NH41_b",
    "MAMNO32":   "NO32_b",
    "MAMNH42":   "NH42_b",
    "MAMNH43":   "NH43_b",
    "MAMSO42":   "SO42_b",
    "MAMSO43":   "SO43_b",
    "MAMCBSO41": "CBSO41_b",
}
MAM_AFTER_KEYS = {
    "HNO3":      "HNO3_a",
    "NH3":        "NH3_a",
    "MAMSO41":   "SO41_a",
    "MAMNO31":   "NO31_a",
    "MAMNH41":   "NH41_a",
    "MAMNO32":   "NO32_a",
    "MAMNH42":   "NH42_a",
    "MAMNH43":   "NH43_a",
    "MAMSO42":   "SO42_a",
    "MAMSO43":   "SO43_a",
    "MAMCBSO41": "CBSO41_a",
    "Nu1":       "Nu1_a",
    "Nu2":       "Nu2_a",
}

def parse_mam(path):
    steps = []
    cur = {}
    in_before = False
    in_after  = False
    step_count = 0

    with open(path) as f:
        for line in f:
            if "MAM_DEBUG_DRIV BEFORE:" not in line and "MAM_DEBUG_DRIV AFTER:" not in line:
                continue

            if "MAM_DEBUG_DRIV BEFORE:" in line:
                if in_after and cur:          # flush previous step
                    steps.append(cur)
                    cur = {}
                in_before = True
                in_after  = False
                # Parse key
                part = line.split("MAM_DEBUG_DRIV BEFORE:")[-1].strip()
                for k, v in MAM_BEFORE_KEYS.items():
                    if part.startswith(k):
                        val = grab(line, k)
                        if val is not None:
                            cur[v] = val
                        break

            elif "MAM_DEBUG_DRIV AFTER:" in line:
                in_before = False
                in_after  = True
                part = line.split("MAM_DEBUG_DRIV AFTER:")[-1].strip()
                for k, v in MAM_AFTER_KEYS.items():
                    if part.startswith(k):
                        val = grab(line, k)
                        if val is not None:
                            cur[v] = val
                        break

    # flush last
    if cur:
        steps.append(cur)
    return steps

# ──────────────────────────────────────────────────────────────────────────────
# STD log parser (HETP)
# ──────────────────────────────────────────────────────────────────────────────
STD_BEFORE_KEYS = {"RH": "RH", "HNO3": "HNO3_b", "NH3": "NH3_b",
                   "SO4": "SO4_b", "NIT": "NIT_b", "NH4": "NH4_b"}
STD_AFTER_KEYS  = {"HNO3": "HNO3_a", "NH3": "NH3_a",
                   "SO4": "SO4_a", "NIT": "NIT_a", "NH4": "NH4_a"}

def parse_std(path):
    steps = []
    cur = {}
    in_after = False

    with open(path) as f:
        for line in f:
            if "MAM_DEBUG_HETP BEFORE:" in line:
                if in_after and cur:
                    steps.append(cur)
                    cur = {}
                in_after = False
                part = line.split("MAM_DEBUG_HETP BEFORE:")[-1].strip()
                for k, v in STD_BEFORE_KEYS.items():
                    if part.startswith(k):
                        val = grab(line, k)
                        if val is not None:
                            cur[v] = val
                        break

            elif "MAM_DEBUG_HETP AFTER:" in line:
                in_after = True
                part = line.split("MAM_DEBUG_HETP AFTER:")[-1].strip()
                for k, v in STD_AFTER_KEYS.items():
                    if part.startswith(k):
                        val = grab(line, k)
                        if val is not None:
                            cur[v] = val
                        break

    if cur:
        steps.append(cur)
    return steps

# ──────────────────────────────────────────────────────────────────────────────
# Main analysis
# ──────────────────────────────────────────────────────────────────────────────

def fnh4_free1(d, which="b"):
    """Free NH4 in mode 1 = NH41 - 2*SO41*(18/96)  [positive = surplus]"""
    nh = d.get(f"NH41_{which}", math.nan)
    so = d.get(f"SO41_{which}", math.nan)
    return nh - 2.0*so*(18.0/96.0)

def fnh4_free2(d, which="b"):
    """Free NH4 in mode 2 = NH42 - SO42*(18/96) (Aitken mode, mostly (NH4)2SO4 too)"""
    nh = d.get(f"NH42_{which}", math.nan)
    so = d.get(f"SO42_{which}", math.nan)
    return nh - 2.0*so*(18.0/96.0)

def safe(d, k):
    return d.get(k, math.nan)

mam = parse_mam(MAM_LOG)
std = parse_std(STD_LOG)

print(f"MAM steps parsed: {len(mam)}")
print(f"STD steps parsed: {len(std)}")
print()

# ── Table 1: Multi-mode SO4 & NH4 balance (BEFORE values each step) ──────────
print("=" * 120)
print(f"{'':>4s} {'RH':>5s} | {'SO41_b':>10s} {'SO42_b':>10s} {'SO43_b':>10s} | "
      f"{'NH41_b':>10s} {'NH42_b':>10s} {'NH43_b':>10s} | "
      f"{'fNH4_1':>10s} {'fNH4_2':>10s} | "
      f"{'NO31_b':>10s} {'NO32_b':>10s}")
print("-" * 120)
hdr = ("Step |  RH   | SO4_m1(kg) SO4_m2(kg) SO4_m3(kg) | "
       "NH4_m1(kg) NH4_m2(kg) NH4_m3(kg) | freeNH4_1  freeNH4_2 | NO3_m1(kg) NO3_m2(kg)")
print(hdr)
print("-" * 120)

for i, d in enumerate(mam, 1):
    fn1 = fnh4_free1(d, "b")
    fn2 = fnh4_free2(d, "b")
    flag = ""
    if not math.isnan(fn1) and fn1 < 0:
        flag = " *** DEFICIT"
    print(f"{i:4d} {safe(d,'RH'):5.3f} | "
          f"{safe(d,'SO41_b'):10.1f} {safe(d,'SO42_b'):10.1f} {safe(d,'SO43_b'):10.1f} | "
          f"{safe(d,'NH41_b'):10.1f} {safe(d,'NH42_b'):10.1f} {safe(d,'NH43_b'):10.1f} | "
          f"{fn1:10.1f} {fn2:10.1f} | "
          f"{safe(d,'NO31_b'):10.1f} {safe(d,'NO32_b'):10.1f}{flag}")

print()
print("=" * 120)
print("MAM AFTER values: SO4, NO3 (total modes 1+2+3) and NH4 AFTER each amicphys call")
print("-" * 120)
hdr2 = ("Step |  RH   | SO41_a     SO42_a     SO43_a  | "
        "NO31_a     NO32_a     NO3tot | "
        "NH41_a     NH42_a     NH43_a | "
        "freeNH4_1a | CBSO41_a  |  Nu1_a       Nu2_a")
print(hdr2)
print("-" * 120)
for i, d in enumerate(mam, 1):
    no3t = safe(d,'NO31_a') + safe(d,'NO32_a') if not (math.isnan(safe(d,'NO31_a')) or math.isnan(safe(d,'NO32_a'))) else math.nan
    fn1a = fnh4_free1(d, "a")
    print(f"{i:4d} {safe(d,'RH'):5.3f} | "
          f"{safe(d,'SO41_a'):10.1f} {safe(d,'SO42_a'):10.1f} {safe(d,'SO43_a'):10.1f} | "
          f"{safe(d,'NO31_a'):10.1f} {safe(d,'NO32_a'):10.1f} {no3t:10.1f} | "
          f"{safe(d,'NH41_a'):10.1f} {safe(d,'NH42_a'):10.1f} {safe(d,'NH43_a'):10.1f} | "
          f"{fn1a:10.1f} | "
          f"{safe(d,'CBSO41_a'):10.1f} | "
          f"{safe(d,'Nu1_a'):12.3e} {safe(d,'Nu2_a'):12.3e}")

print()
print("=" * 120)
print("STD HETP (ISORROPIA) AFTER values")
print("-" * 120)
print(f"{'Step':>4s} {'RH':>5s} | {'SO4_a':>10s} {'NIT_a':>10s} {'NH4_a':>10s} | "
      f"{'HNO3_a':>10s} {'NH3_a':>10s} | {'NIT_b':>10s}")
print("-" * 120)
for i, d in enumerate(std, 1):
    print(f"{i:4d} {safe(d,'RH'):5.3f} | "
          f"{safe(d,'SO4_a'):10.1f} {safe(d,'NIT_a'):10.1f} {safe(d,'NH4_a'):10.1f} | "
          f"{safe(d,'HNO3_a'):10.1f} {safe(d,'NH3_a'):10.1f} | "
          f"{safe(d,'NIT_b'):10.1f}")

print()
print("=" * 120)
print("MAM vs STD NO3 partitioning fraction  [after each chemistry call]")
print("  f_MAM = (NO31+NO32) / (NO31+NO32+HNO3)   [gas+aero_MAM total]")
print("  f_STD = NIT / (NIT + HNO3)                [gas+aero_STD total]")
print("-" * 120)
print(f"{'Step':>4s} {'RH':>5s} | {'MAM NO31':>10s} {'MAM NO32':>10s} {'MAM HNO3':>10s} {'f_MAM':>8s} | "
      f"{'STD NIT':>10s} {'STD HNO3':>10s} {'f_STD':>8s} | ratio f_MAM/f_STD")
print("-" * 120)
nsteps = min(len(mam), len(std))
for i in range(nsteps):
    dm, ds = mam[i], std[i]
    n31 = safe(dm,'NO31_a')
    n32 = safe(dm,'NO32_a')
    hno3m = safe(dm,'HNO3_a')
    tot_m = n31 + n32 + hno3m
    f_mam = (n31+n32)/tot_m if tot_m > 0 else math.nan

    nit_s = safe(ds,'NIT_a')
    hno3s = safe(ds,'HNO3_a')
    tot_s = nit_s + hno3s
    f_std = nit_s/tot_s if tot_s > 0 else math.nan

    ratio = f_mam/f_std if (not math.isnan(f_mam) and not math.isnan(f_std) and f_std>0) else math.nan
    rh = safe(dm,'RH')
    print(f"{i+1:4d} {rh:5.3f} | "
          f"{n31:10.1f} {n32:10.1f} {hno3m:10.1f} {f_mam:8.4f} | "
          f"{nit_s:10.1f} {hno3s:10.1f} {f_std:8.4f} | {ratio:8.4f}")

print()
# Focus on the transition around step 47
print("=" * 120)
print("STEP-47 ZOOM: SO4 change between AFTER(step N-1) and BEFORE(step N)")
print("  [AFTER - BEFORE delta reveals external SO4 injection between calls]")
print("-" * 120)
print(f"{'Step':>4s} | {'SO41_b':>10s} {'SO41_a-prev':>12s} {'dSO41':>10s} | "
      f"{'SO42_b':>10s} {'SO42_a-prev':>12s} {'dSO42':>10s} | "
      f"{'NH41_b':>10s} {'NH41_a-prev':>12s} {'dNH41':>10s}")
print("-" * 120)
for i in range(1, min(len(mam), 60)):
    prev = mam[i-1]
    curr = mam[i]
    d1 = safe(curr,'SO41_b') - safe(prev,'SO41_a')
    d2 = safe(curr,'SO42_b') - safe(prev,'SO42_a')
    dn = safe(curr,'NH41_b') - safe(prev,'NH41_a')
    flag = " <-- JUMP" if abs(d1)>100 or abs(d2)>100 else ""
    print(f"{i+1:4d} | "
          f"{safe(curr,'SO41_b'):10.1f} {safe(prev,'SO41_a'):12.1f} {d1:10.1f} | "
          f"{safe(curr,'SO42_b'):10.1f} {safe(prev,'SO42_a'):12.1f} {d2:10.1f} | "
          f"{safe(curr,'NH41_b'):10.1f} {safe(prev,'NH41_a'):12.1f} {dn:10.1f}{flag}")
