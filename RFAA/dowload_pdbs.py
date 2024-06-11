import requests

def download_pdb(pdb_id, output_dir="pdb_files"):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if response.status_code == 200:
        with open(f"{output_dir}/{pdb_id}.pdb", 'w') as file:
            file.write(response.text)
        print(f"Downloaded: {pdb_id}")
    else:
        print(f"Failed to download: {pdb_id}")

def main(pdb_list_file, output_dir="pdb_files"):
    # Create the output directory if it doesn't exist
    import os
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with open(pdb_list_file, 'r') as file:
        pdb_ids = file.readlines()
    
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.strip()
        if pdb_id:
            download_pdb(pdb_id, output_dir)

if __name__ == "__main__":
    pdb_list_file = "pdb_names.txt"
    main(pdb_list_file)
