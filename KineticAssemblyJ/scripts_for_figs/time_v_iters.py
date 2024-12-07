import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

# Initialize a nested dictionary to hold the data
# Format: data[integrator][nmer_size] = peak_memory_usage
for_iters = []
for_time = []
rev_iters = []
rev_time = []
py_iters = [1,10,100,1000,2000,3000]
py_time = [0.38211774826049805,4.842249631881714,51.591909408569336,480.18445110321045,1033.9070744514465,1591.7473483085632]
directory = "/home/stakesh1/scr4_mjohn218/sho/KineticAssembly/KineticAssemblyJ/bench_iters"
for filename in os.listdir(directory):
    if filename.endswith('.out'):
        
        filepath = os.path.join(directory, filename)
        with open(filepath, 'r') as file:    
            lines = file.readlines()
            if len(lines)<8:
                continue
            time = 0
            if lines[-2].strip().split()[1] == "ms":
                time = float(lines[-2].strip().split()[0])/1000
            else:
                time = float(lines[-2].strip().split()[0])
            
            iters = int(lines[1].strip().split()[1])

            if lines[0].strip() == 'AutoForwardDiff()':
                for_iters.append(iters)
                for_time.append(time)
            else:
                rev_iters.append(iters)
                rev_time.append(time)



for_df = pd.DataFrame({'Iterations': for_iters, 'Time (s)': for_time, 'Method': 'ForwardDiff'})
rev_df = pd.DataFrame({'Iterations': rev_iters, 'Time (s)': rev_time, 'Method': 'ReverseDiff'})
python = pd.DataFrame({'Iterations': py_iters, 'Time (s)': py_time, 'Method': 'Python'})
df = for_df.append(rev_df)
df = df.append(python)
print(df)

plt.figure(figsize=(8, 6))
plt.style.use('seaborn-white')
plt.xscale('log')
plt.xlabel('Iterations',fontsize=20)
plt.ylabel('Time (s)',fontsize=20)

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
#plot = sns.lineplot(data=df, x="Iterations", y="Time (s)",hue='Method')
plot = sns.lineplot(data=for_df, x="Iterations", y="Time (s)",label='ForwardDiff')
plt.setp(plot.get_legend().get_texts(), fontsize='15')

#plot.get_figure().savefig('time_v_iters.png', format="png")
plot.get_figure().savefig('time_v_iters_forw.png', format="png")