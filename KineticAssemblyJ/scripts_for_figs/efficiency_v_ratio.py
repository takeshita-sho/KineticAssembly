import re
import matplotlib.pyplot as plt

# Initialize lists to hold the values
ktri_kdim_values = []
efficiency_values = []

# Open and read the file
with open('optim_job_18479658.out', 'r') as file:
    lines = file.readlines()

# Flag to start parsing after the header line
parsing = False

for line in lines:
    line = line.strip()
    
    # Start parsing after the header
    if line == 'Iteration  Ktri/Kdim  Efficiency':
        parsing = True
        continue
    if parsing:
        # Match lines with data
        match = re.match(r'(\d+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)', line)
        if match:
            ktri_kdim = float(match.group(2))
            efficiency = float(match.group(3))
            ktri_kdim_values.append(ktri_kdim)
            efficiency_values.append(efficiency)
        else:
            # Stop parsing if the line doesn't match
            parsing = False

# Plot the data
plt.style.use('seaborn-white')
plt.figure(figsize=(8, 6))
plt.scatter(ktri_kdim_values, efficiency_values, marker='o',s=100)
plt.xlabel('Ktri/Kdim',size=20)
plt.ylabel('Efficiency',size=20)
plt.title('Efficiency vs Ktri/Kdim for a Heptamer',size=20)
plt.xticks(size=20)
plt.yticks(size=20)
#plt.xlim(3, 5.5) 
plt.tight_layout()
plt.savefig('efficiency_vs_ktri_kdim_7mer.png')