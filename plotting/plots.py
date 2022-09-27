import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from matplotlib import rc
rc("text", usetex=True)


def load_data(filename):
    data = pd.read_csv(filename, delimiter=" ", header=0, skiprows=1)
    data['T'] = 1/data['beta']
    f = open(filename)
    n = f.readline()
    n = int(n[2:])

    return n, data


def average_over_cycles(data):
    data['Energy_squared'] = data['Energy'] * data['Energy']
    avg = data.groupby(by=['beta'], as_index=False).mean()
    avg['T'] = 1/avg['beta']
    betas = 1 / avg.index
    return avg

def std_over_cycles(data):
    std = data.groupby(by=['beta'], as_index=False).std()
    std['T'] = 1 / std['beta']
    return std


def plot(data):
    #n, data = load_data(filename)
    #calcP2s(data)


    avg = average_over_cycles(data)


    # plt.xlabel('T')
    # plt.ylabel('P2, P2m')
    # plt.plot(avg['T'], avg['P2'])
    # plt.plot(avg['T'], avg['P2m'])
    # plt.show()
    plt.subplot(3,1,1)
    plt.xlabel('T')
    plt.ylabel('E')
    plt.plot(avg['T'], avg['Energy'])
    plt.scatter(data['T'], data['Energy'], s=2)
    plt.plot(avg['T'], avg['Energy'])
    plt.subplot(3,1,2)
    plt.xlabel('T')
    plt.ylabel('Cv')
    plt.plot(avg['T'], n**3*(-avg['Energy']**2+avg['Energy_squared'])*avg['beta']**2)
    plt.subplot(3, 1, 3)
    plt.xlabel('T')
    plt.ylabel('Cv')
    plt.plot(avg['T'], n ** 3 * (-avg['Energy'] ** 2 + avg['Energy_squared']) * avg['beta'] ** 2)
    #plt.plot(avg['T'], np.gradient(avg['Energy'], avg['T']))



def calcP2(data):
    split = [y for x, y in data.groupby(['beta'], as_index=False)]

    betas = []
    P2 = np.zeros((3,3))
    P2s = []
    P2ms = []
    Nsample = 500
    for group in split:
        
        betas.append(group['beta'][group.index.min()])
        randints = random.sample(set(group.index.values), Nsample)
        p2temp = []
        p2mtemp = []
        for i in range(Nsample):
            P2[0, 0] = group['P2_00'][randints[i]]
            P2[0, 1] = group['P2_01'][randints[i]]
            P2[0, 2] = group['P2_02'][randints[i]]

            P2[1, 0] = group['P2_10'][randints[i]]
            P2[1, 1] = group['P2_11'][randints[i]]
            P2[1, 2] = group['P2_12'][randints[i]]

            P2[2, 0] = group['P2_20'][randints[i]]
            P2[2, 1] = group['P2_21'][randints[i]]
            P2[2, 2] = group['P2_22'][randints[i]]

            eigvals, eigvecs = np.linalg.eig(P2)
            eigvals = np.sort(eigvals)
            p2temp.append(eigvals[-1])
            p2mtemp.append(eigvals[1])
        P2s.append(np.average(p2temp))
        P2ms.append(np.average(p2mtemp))
    # plt.plot(1/np.array(betas), P2s)
    # plt.plot(1/np.array(betas), P2ms)
    # plt.show()
    return 1/np.array(betas), P2s, P2ms





files =  [#"test.txt"
         "cooldown_zoom_cluster.txt",
         "heatup_zoom_cluster.txt",
         #"heatup_zoom.txt"
         # "00heatup.txt",
         # "01heatup.txt",
         # "02heatup.txt",
         # "03heatup.txt",
         # "04heatup.txt",
         # "05heatup.txt",
         # "06heatup.txt",
         # "07heatup.txt",
         # "08heatup.txt",
         # "09heatup.txt",
         # "00cooldown.txt",
         # "01cooldown.txt",
         # "02cooldown.txt",
         # "03cooldown.txt",
         # "04cooldown.txt",
         # "05cooldown.txt",
         # "06cooldown.txt",
         # "07cooldown.txt",
         # "08cooldown.txt",
         # "09cooldown.txt"
         ]

folder = "../output/"

dataset =[]

labels=["$r=0.0$ heatup",
        "$r=0.1$ heatup",
        "$r=0.2$ heatup",
        "$r=0.3$ heatup",
        "$r=0.4$ heatup",
        "$r=0.5$ heatup",
        "$r=0.6$ heatup",
        "$r=0.7$ heatup",
        "$r=0.8$ heatup",
        "$r=0.9$ heatup",
        "$r=0.0$ cooldown",
        "$r=0.1$ cooldown",
        "$r=0.2$ cooldown",
        "$r=0.3$ cooldown",
        "$r=0.4$ cooldown"
        "$r=0.0$ cooldown",
        "$r=0.1$ cooldown",
        "$r=0.2$ cooldown",
        "$r=0.3$ cooldown",
        "$r=0.4$ cooldown"
        ]
#join two heatup and cooldown data

for file in files:
    n, data = load_data(folder + file)
    dataset.append(data)

dataset.append(pd.concat([dataset[0], dataset[1]], ignore_index=True, sort=False))

print(dataset)


avgs = []
stds = []
P2s = []
P2ms = []
Ts = []
#calculate averages, standard deviations and P2, P2m
for data in dataset:
    avgs.append(average_over_cycles(data))
    stds.append(std_over_cycles(data))
    T, P2, P2m = calcP2(data)
    P2s.append(P2)
    P2ms.append(P2m)
    Ts.append(T)




#plot
for i in range(len(avgs)):
    plt.errorbar(avgs[i]['T'], avgs[i]['Energy'], stds[i]['Energy'], label = labels[i])
    #plt.errorbar(avgs[i]['T'], avgs[i]['polar_order'], stds[i]['polar_order'], label=labels[i])

plt.xlabel('$T$')
plt.ylabel('$E$')
# plt.legend()
plt.show()

for i in range(len(avgs)):
    plt.scatter(avgs[i]['T'], n ** 3 * (-avgs[i]['Energy'] ** 2 + avgs[i]['Energy_squared']) * avgs[i]['beta'] ** 2, label = labels[i])

plt.xlabel('$T$')
plt.ylabel('$C_v$')
# plt.legend()
plt.show()

for i in range(len(P2s)):
    plt.plot(np.array(Ts[i]), P2s[i], label = labels[i])
    plt.plot(np.array(Ts[i]), P2ms[i], label=labels[i])
    # plt.plot(np.array(Ts[i]), P2ms[i], label=labels[i])

plt.xlabel('$T$')
plt.ylabel('$P_2$, $P_{2m}$')
# plt.legend()
plt.show()






