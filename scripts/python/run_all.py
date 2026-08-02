"""
run_all.py -- regenerate every figure and table of
  Monjo (2025), ApJ 981, 195, from data/ using the Python pipeline.
Usage:  python run_all.py
"""
import runpy
import os

HERE = os.path.dirname(os.path.abspath(__file__))
for script in ["fig0_four_clusters.py", "fig1_complete.py", "fig2_table2.py"]:
    print(f"\n===== running {script} =====")
    runpy.run_path(os.path.join(HERE, script), run_name="__main__")
print("\n== DONE. Outputs (*_python.*) in figures/. ==")
