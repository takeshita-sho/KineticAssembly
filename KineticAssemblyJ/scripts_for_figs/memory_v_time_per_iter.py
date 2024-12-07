import re
import matplotlib.pyplot as plt

with open("/home/stakesh1/scr4_mjohn218/sho/KineticAssembly/KineticAssemblyJ/benchmarkiters17940718.out", 'r') as file:
    data={}
    current_iter = None
    for line in file:
        
        iter_match = re.search(r'Iterations\s*(\d+)', line)
        if iter_match:
            current_iter = int(iter_match.group(1))
            if current_iter not in data:
                data[current_iter] = []
        mem_match = re.search(r'Peak memory usage so far:\s*([\d:]+)', line)
        if mem_match:
            current_time = float(mem_match.group(1))
            data[current_iter].append(current_time)

plt.figure(figsize=(10, 6))
#for iterations, values in sorted(data.items()):
for iterations, values in data.items():
    plt.plot(values, label=f'Iteration {iterations}')
plt.xlabel('Time (Loss Calls)')
plt.ylabel('Memory Usage (MB)')
plt.title('Memory Usage vs Time')
plt.legend()
plt.tight_layout()
plt.savefig('mem_v_iters.png')
