import requests
import json

def fetch_antibody(drug_name):
    url = f"https://absd.pasteur.cloud/api/antibodies?search={drug_name}"
    
    print(f"[fetching] {drug_name}")
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"[error] Could not fetch {drug_name}")
        return
    
    data = response.json()
    
    if not data or len(data) == 0:
        print(f"[error] No results for {drug_name}")
        return
    
    print(f"[found] {len(data)} results\n")
    
    for i, antibody in enumerate(data[:3]):
        print(f"--- Result {i+1} ---")
        print(f"Name: {antibody.get('name', 'N/A')}")
        print(f"Heavy Chain: {antibody.get('heavy_chain', 'N/A')[:100]}...")
        print(f"Light Chain: {antibody.get('light_chain', 'N/A')[:100]}...")
        print()

if __name__ == "__main__":
    print("Enter antibody name")
    drug = input("> ").strip()
    fetch_antibody(drug)