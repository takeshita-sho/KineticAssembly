import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

# Initialize a nested dictionary to hold the data
# Format: data[integrator][nmer_size] = peak_memory_usage
for_iters = []
for_mem = []
rev_iters = []
rev_mem = []
py_iters = [1,10,100,1000,2000]
py_mem = [353.18359375,489.26171875,1827.5078125,15092.515625,30093.13671875]
for i in range(len(py_mem)):
    py_mem[i] = py_mem[i]*1.04858
directory = "/home/stakesh1/scr4_mjohn218/sho/KineticAssembly/KineticAssemblyJ/bench_iters"
for filename in os.listdir(directory):
    if filename.endswith('.out'):
        
        filepath = os.path.join(directory, filename)
        with open(filepath, 'r') as file:    
            lines = file.readlines()
            if len(lines)<8:
                continue
            mem = float(lines[-3].strip().split()[0])
            
            iters = int(lines[1].strip().split()[1])

            if lines[0].strip() == 'AutoForwardDiff()':
                for_iters.append(iters)
                for_mem.append(mem)
            else:
                rev_iters.append(iters)
                rev_mem.append(mem)



for_df = pd.DataFrame({'Iterations': for_iters, 'Peak Memory (MB)': for_mem, 'Method': 'ForwardDiff'})
rev_df = pd.DataFrame({'Iterations': rev_iters, 'Peak Memory (MB)': rev_mem, 'Method': 'ReverseDiff'})
python = pd.DataFrame({'Iterations': py_iters, 'Peak Memory (MB)': py_mem, 'Method': 'Python'})
df = pd.concat([for_df, rev_df, python], ignore_index=True)
average_forward_mem = for_df.groupby('Iterations')['Peak Memory (MB)'].mean().reset_index()
average_reverse_mem = rev_df.groupby('Iterations')['Peak Memory (MB)'].mean().reset_index()
average_python_mem = python.groupby('Iterations')['Peak Memory (MB)'].mean().reset_index()

merged_df = pd.merge(
    average_forward_mem,
    average_reverse_mem,
    on='Iterations',
    how='inner'
)

# Merge the result with Python
merged_df = pd.merge(
    merged_df,
    average_python_mem,
    on='Iterations',
    how='inner'
)
merged_df.rename(columns={
    'Peak Memory (MB)_x': 'ForwardDiff',
    'Peak Memory (MB)_y': 'ReverseDiff',
    'Peak Memory (MB)': 'Python'
}, inplace=True)
print(merged_df['ForwardDiff']/merged_df['Python'])
#print(merged_df['Peak Memory (MB)'])

# Pivot the DataFrame to have Methods as separate columns
pivot_df = df.pivot(index='Iterations', columns='Method', values='Peak Memory (MB)')

# Drop any rows with missing values
pivot_df = pivot_df.dropna()

# Calculate the ratios
pivot_df['Reverse/Forward'] = pivot_df['ReverseDiff'] / pivot_df['ForwardDiff']
pivot_df['Python/Forward'] = pivot_df['Python'] / pivot_df['ForwardDiff']

# Reset index for melting
pivot_df_reset = pivot_df.reset_index()

# Melt the DataFrame to long format
ratio_df = pivot_df_reset.melt(
    id_vars='Iterations',
    value_vars=['Reverse/Forward', 'Python/Forward'],
    var_name='Ratio_Type',
    value_name='Ratio'
)
plot = sns.lineplot(data=ratio_df, x="Iterations", y="Ratio",hue='Ratio_Type')
plt.figure(figsize=(8, 6))
plt.style.use('seaborn-white')
plt.xscale('log')
plt.xlabel('Iterations',fontsize=20)
plt.ylabel('Peak Memory (MB)',fontsize=20)

plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.tight_layout()
plt.setp(plot.get_legend().get_texts(), fontsize='15')
plot.get_figure().savefig('mem_v_iters.png', format="png")        
            
