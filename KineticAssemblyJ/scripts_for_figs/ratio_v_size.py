import matplotlib.pyplot as plt
import numpy as np


k_tri_fully_connected = []
# Data for the plot
for i in range(1,5):
    with open(f'optim_job_18357681_{i}.out', 'r') as file:
        lines = file.readlines()
        max_eff = 0
        max_ratio = 0
        parsing = False
        for l in lines:
            line = l.strip()
            if line == 'Iteration  Ktri/Kdim  Efficiency':
                parsing = True
                continue
            if parsing:
                # Match lines with data
                args = line.split()
                
                eff = float(args[2])
                ratio = float(args[1])
                if eff > max_eff:
                    max_eff = eff
                    max_ratio = ratio
        k_tri_fully_connected.append(max_ratio)  
with open('optim_job_18479658.out', 'r') as file:
        lines = file.readlines()
        max_eff = 0
        max_ratio = 0
        parsing = False
        for l in lines:
            line = l.strip()
            if line == 'Iteration  Ktri/Kdim  Efficiency':
                parsing = True
                continue
            if parsing:
                # Match lines with data
                args = line.split()
                
                eff = float(args[2])
                ratio = float(args[1])
                if eff > max_eff:
                    max_eff = eff
                    max_ratio = ratio
        k_tri_fully_connected.append(max_ratio)  

assembly_size = [3, 4, 5, 6,7]

# Plotting
plt.figure(figsize=(8, 8))
plt.plot(assembly_size, k_tri_fully_connected, 
         color='forestgreen', marker='o',  label='Fully Connected',lw=7,markersize=20)


# Logarithmic scale for the y-axis
plt.yscale('log')

# Labels and title
plt.xlabel('Assembly Size', fontsize=40)
plt.ylabel(r'$k_{tri} / k_{dim}$', fontsize=40)
#plt.title('Comparison of Network Topologies', fontsize=14)
print(k_tri_fully_connected)
# Legend
#plt.legend(loc='upper left', fontsize=20,frameon=False)
plt.ylim(top=78)
#plt.xlim(right=7.3)
# Custom ticks
plt.xticks([3,4,5,6,7], fontsize=40)
plt.yticks(fontsize=40)
plt.tick_params(axis='both', which='both',width=3,length=10)

# Grid and light background
#plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)


# Display the plot
plt.tight_layout()
plt.savefig('ratio_v_size.png', format="png")