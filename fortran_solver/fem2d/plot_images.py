import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Configuration
input_dir = '.'          # Directory containing the P_*.txt files
output_dir = 'images'    # Directory to save plots
file_prefix = 'P_'       # Prefix of files to process
file_suffix = '.txt'     # Suffix of files to process

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Set better plot style
rcParams['figure.figsize'] = (8, 8)
rcParams['axes.grid'] = True
rcParams['grid.alpha'] = 0.3
rcParams['axes.axisbelow'] = True

# Get all matching files in the directory
file_list = [f for f in os.listdir(input_dir) 
             if f.startswith(file_prefix) and f.endswith(file_suffix)]
file_list.sort()  # Sort files numerically

print(f"Found {len(file_list)} files to process")

# Process each file
for filename in file_list:
    # Extract the base name for plot title and output filename
    basename = os.path.splitext(filename)[0]
    
    # Load data (assuming two rows: x and y coordinates)
    try:
        data = np.loadtxt(os.path.join(input_dir, filename))
        
        # Check if we have exactly two rows (x and y)
        if data.shape[0] != 2:
            print(f"Skipping {filename}: expected 2 rows, found {data.shape[0]}")
            continue
            
        x_coords = data[0, :]
        y_coords = data[1, :]
        
        # Create plot
        plt.figure()
        plt.scatter(x_coords, y_coords, s=10, alpha=0.7, c='blue')
        plt.title(f'Point Cloud: {basename}')
        plt.xlabel('X coordinate')
        plt.ylabel('Y coordinate')
        plt.axis('equal')  # Keep aspect ratio equal
        
        # Save plot
        output_path = os.path.join(output_dir, f'{basename}.png')
        plt.xlim([-0.001,0.001])
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Processed {filename} -> {output_path}")
        
    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")

print("Processing complete!")
