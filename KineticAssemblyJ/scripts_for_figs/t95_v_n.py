import matplotlib.pyplot as plt
import numpy as np

t95 = []
# Data for the plot
for i in range(1,6):
    with open(f't95_18390519_{i}.out', 'r') as file:
        lines = file.readlines()
        min_t95 = float('inf')
        parsing = False
        for l in lines:
            line = l.strip()
            if line == 'Iteration  t95':
                parsing = True
                continue
            if parsing:
                # Match lines with data
                args = line.split()
                
                curr_t95 = float(args[1])
                
                if 0 < curr_t95 < min_t95:
                    min_t95 = curr_t95
        t95.append(min_t95)  

assembly_size = [3, 4, 5, 6,7]

# Plotting
plt.figure(figsize=(8, 8))
plt.plot(assembly_size, t95, 
         color='forestgreen', marker='o',  label='Fully Connected',lw=7,markersize=20)


# Logarithmic scale for the y-axis
plt.yscale('log')

# Labels and title
plt.xlabel('Assembly Size', fontsize=40)
plt.ylabel(r'$\tau^{*}_{95}$', fontsize=40)
#plt.title('Comparison of Network Topologies', fontsize=14)

# Legend
#plt.legend(loc='upper left', fontsize=30,frameon=False)
plt.ylim(top=130,bottom=50)
# Custom ticks
plt.xlim(right=7.1)
plt.xticks([3,4,5,6,7], fontsize=40)
plt.yticks(fontsize=40)
plt.tick_params(axis='both', which='both',width=3,length=10)

# Grid and light background
#plt.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)


# Display the plot
plt.tight_layout()
plt.savefig('t95_v_size.png', format="png")