from spmUtils import *


folder_path = "/home/captainbroccoli/Documents/2022-06-30/"
prefix = "20220630-152823_Cu(111)--AFM_NonContact_QPlus_AtomManipulation--"
sufix = ".Df(Z)_mtrx"

path_on = f"{folder_path}{prefix}130_1{sufix}"
path_off = f"{folder_path}{prefix}134_1{sufix}"

mtrx_data = access2thematrix.MtrxData()
data_file = f'{path_off}'
traces, message = mtrx_data.open(data_file)
curve_trace, message = mtrx_data.select_curve(traces[0])

print(curve_trace.referenced_by["Data File Name"])

x = curve_trace.data[0]
y = curve_trace.data[1]

fig = plt.figure()
sns.lineplot(x,y)
plt.show()
