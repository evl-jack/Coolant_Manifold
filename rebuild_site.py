from pathlib import Path
import subprocess,sys
r=Path(__file__).resolve().parent
for script in ('export_results.py','build_viewer.py','build_site.py'):
    subprocess.run([sys.executable,str(r/script)],cwd=r,check=True)
