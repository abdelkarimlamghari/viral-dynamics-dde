import subprocess
import sys

scripts = [
    'scripts/run_main.py',
    'scripts/run_hopf.py',
    'scripts/run_sensitivity.py',
    'scripts/run_clinical.py'
]

for scr in scripts:
    print(f"\n>>> Running {scr}")
    ret = subprocess.run([sys.executable, scr])
    if ret.returncode != 0:
        print(f"ERROR in {scr}")
        sys.exit(1)

print("\n✅ All figures generated successfully in figures/")
