import matplotlib.pyplot as plt
filename = 'turbo_data.txt'
value=[]
with open(filename, 'r') as f:
    lines = f.readlines()
    for line in lines:
        value.append( [float(s) for s in line.split()])

plt.axes(yscale = "log")
plt.title("turbo decoder")
for i in range(len(value)-1):
    plt.plot(value[0],value[i+1],marker='o',linestyle='--')
plt.grid(True,linestyle = "--",color = 'gray' ,linewidth = '0.5',axis='both')

plt.legend(['L100_I10', 'L1024_I10', 'L1024_I18'])
plt.xlabel('SNR (dB)')
plt.ylabel('BER')

plt.savefig('performance.png')
plt.show()

