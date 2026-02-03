import os
import time
from pathlib import Path

# Target file path
FILE_PATH = Path(r"C:\Users\bunsr\rosalind-bioinformatics\multispecific_antibodies\Whole_mAb\mAb_3D_Physics_Features.csv")

print(f"Monitoring for: {FILE_PATH.name}...")

while True:
    if FILE_PATH.exists():
        # Check if the file size is stable (ensures writing is finished)
        size_1 = FILE_PATH.stat().st_size
        time.sleep(2)
        size_2 = FILE_PATH.stat().st_size
        
        if size_1 == size_2 and size_1 > 0:
            print(f"\n[!!!] FILE CREATED SUCCESSFULLY: {FILE_PATH}")
            print(f"Final Size: {size_1 / 1024:.2f} KB")
            break
    
    time.sleep(5)  # Check every 5 seconds to save CPU cycles