import pandas as pd
import os
import lime
# Folder containing the .txt files
folder_path = "/home/vital/Dropbox/Astrophysics/Data/CAPERS/sample/CAPERS_EGS_V0.2.2/bands"

# Target index and column
target_index = 'He1_10832A'
target_column = 'wavelength'
new_value = 10832.0570

for filename in os.listdir(folder_path):
    if filename.endswith(".txt"):
        file_path = os.path.join(folder_path, filename)

        try:

            df = lime.load_frame(file_path)  # Adjust sep if needed

            if target_index in df.index and target_column in df.columns:
                df.at[target_index, target_column] = new_value
                lime.save_frame(file_path, df)
                print(f"Updated {filename}")

        except Exception as e:
            print(f"Could not process {filename}: {e}")