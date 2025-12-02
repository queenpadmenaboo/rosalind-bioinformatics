# Update readme_count.py

import os
import re
from datetime import datetime

def update_readme_count():
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Change to multispecific_antibodies directory
    os.chdir(script_dir)
    
    # Debug: show current directory
    print(f"Current directory: {os.getcwd()}")
    print(f"Files found: {os.listdir('.')}")
    
    # Count .py files (excluding sabdabconverter.py and this script)
    py_files = [f for f in os.listdir('.') 
                if f.endswith('.py') 
                and f not in ['sabdabconverter.py', 'readme_count.py', 'selenium_antibody_scraper.py', 'thera_sabdab_scraper.py']]
    
    print(f"Antibody .py files: {py_files}")
    count = len(py_files)
    
    readme_path = 'README.md'
    if os.path.exists(readme_path):
        with open(readme_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Update the count line
        current_date = datetime.now().strftime("%B %Y")
        
        updated_content = re.sub(
            r'- Data last updated:.*\n- Total antibodies:.*',
            f'- Data last updated: {current_date}\n- Total antibodies: {count}',
            content
        )
        
        with open(readme_path, 'w', encoding='utf-8') as f:
            f.write(updated_content)
        
        print(f"README updated: {count} antibodies, {current_date}")
    else:
        print("README.md not found!")

if __name__ == "__main__":
    update_readme_count()