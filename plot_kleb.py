import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# Check if the CSV file name is provided as a command-line argument
if len(sys.argv) < 2:
    print("Usage: python script_name.py file_name.csv")
    sys.exit(1)

file_path = sys.argv[1]  # The CSV file name is the first argument

# Read the CSV file
df = pd.read_csv(file_path, header=None)

# Extract x and y values from the DataFrame
x_values = df.iloc[0].values
y_values = df.iloc[1].values

# Create a bar plot
plt.figure(figsize=(10, 6))
plt.bar(range(len(x_values)), y_values, color='blue')  # Use range(len(x_values)) as x values
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
plt.title('Bar Plot')

# Set tick labels to be the x values
plt.xticks(range(len(x_values)), x_values, rotation=45)  # Rotate x-axis labels for better readability
plt.grid(True)
plt.tight_layout()  # Adjust layout to prevent clipping of x-axis labels
plt.show()


""" # Load the CSV file
df = pd.read_csv(file_path)

# Transposing the data since it's in a single row
data = df.transpose()

# Removing zero values and compute the logarithm of the data
#log_data = np.log2(data+1)
log_data = data
# Plotting the values as a scatter plot with log scale for non-zero values
plt.figure(figsize=(30, 6))
plt.scatter(log_data.index, log_data[0], color='blue')
plt.title('Scatter Plot of the Log of Non-Zero Values')
plt.xlabel('Index')
plt.ylabel('Log2 Value')
plt.grid(True)
plt.show()
 """