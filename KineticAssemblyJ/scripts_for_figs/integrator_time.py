import re
import pandas as pd
import matplotlib.pyplot as plt
import os

# Initialize a nested dictionary to hold the data
# Format: data[integrator][nmer_size] = peak_memory_usage
data = {}
directory = "/home/stakesh1/scr4_mjohn218/sho/KineticAssembly/KineticAssemblyJ/data"
# Open and read the file
for filename in os.listdir(directory):
    
    if filename.endswith('.out'):
        
        filepath = os.path.join(directory, filename)
        with open(filepath, 'r') as file:
            lines = file.readlines()
            
            if lines[2].strip() == 'AutoForwardDiff()':
                # Match nmer size (e.g., "nmer_size = 5")
                nmer_match = re.match(r'(\d+)mer', lines[4].strip())
                current_nmer = int(nmer_match.group(1))
                if current_nmer < 7:
                
                    # Match integrator names (e.g., "IntegratorName{")
                    integrator_match = re.match(r'(\w+)\{', lines[3].strip())
                    current_integrator = integrator_match.group(1)
                    if current_integrator not in data:
                        data[current_integrator] = {}
                    
                    
                        
                    
                    # Match peak memory usage (e.g., "Peak memory usage so far: 2012.78125")
                    mem_match = lines[-2].strip().split()
                    
                    
                    if mem_match[1] == 'ms':
                        memory_usage = float(mem_match[0])/1000
                    else: 
                        memory_usage = float(mem_match[0])
                    
                    # Update the maximum memory usage for this integrator and nmer size
                    #data[current_integrator][current_nmer] = memory_usage
                    if current_nmer not in data[current_integrator]:
                        data[current_integrator][current_nmer] = memory_usage
                    else:
                        data[current_integrator][current_nmer]+= memory_usage

for integrator in data:
    for nmer_size in data[integrator]:
        data[integrator][nmer_size] /= 3

# Convert the data into a pandas DataFrame for easier plotting
df = pd.DataFrame(data).sort_index()
print(df)
# Plotting
plt.figure(figsize=(8, 6))
plt.style.use('seaborn-white')
for integrator in df.columns:
    plt.plot(df.index, df[integrator], marker='o', label=integrator,linewidth=5,markersize=10)

plt.xlabel('Nmer Size',size=20)
plt.ylabel('Time (s)',size=20)
plt.xticks(range(3,7),size=20)
plt.yticks(size=20)
#plt.ylim(0,120)
#plt.xticks(range(3,9))
plt.title('ForwardDiff',size=20)
plt.legend(fontsize=15)
plt.tight_layout()
plt.savefig('integrator_time_for.png')
plt.show()
