import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


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

def calcP2s(data):
    P2s = []
    P2ms = []
    P2 = np.zeros((3,3))
    for i in range(data['P2_00'].size):
        P2[0, 0] = data['P2_00'][i]
        P2[0, 1] = data['P2_01'][i]
        P2[0, 2] = data['P2_02'][i]

        P2[1, 0] = data['P2_10'][i]
        P2[1, 1] = data['P2_11'][i]
        P2[1, 2] = data['P2_12'][i]

        P2[2, 0] = data['P2_20'][i]
        P2[2, 1] = data['P2_21'][i]
        P2[2, 2] = data['P2_22'][i]

        eigvals, eigvecs = np.linalg.eig(P2)
        eigvals = np.sort(eigvals)
        P2s.append(eigvals[-1])
        P2ms.append(eigvals[1])

    data['P2'] = P2s
    data['P2m'] = P2ms

    return data

def plot(filename):
    n, data = load_data(filename)
    #calcP2s(data)

    avg = average_over_cycles(data)

    # plt.xlabel('T')
    # plt.ylabel('P2, P2m')
    # plt.plot(avg['T'], avg['P2'])
    # plt.plot(avg['T'], avg['P2m'])
    # plt.show()
    plt.subplot(2,1,1)
    plt.xlabel('T')
    plt.ylabel('E')
    plt.plot(avg['T'], avg['Energy'])
    plt.scatter(data['T'], data['Energy'], s=2)
    plt.plot(avg['T'], avg['Energy'])
    plt.subplot(2,1,2)
    plt.xlabel('T')
    plt.ylabel('Cv')
    plt.plot(avg['T'], n**3*(-avg['Energy']**2+avg['Energy_squared'])*avg['beta']**2)
    #plt.plot(avg['T'], np.gradient(avg['Energy'], avg['T']))



def plotF(filename):
    n, data = load_data(filename)
    for i in range(50):

       plt.title(data.iloc[i*1000]['T'])
       plt.hist(data[i*1000:(i+1)*1000-1]['Energy'], bins=100)
       plt.show()






#plotF("cooldown_zoom_cluster.txt")

#plot("heatup_zoom.txt")
#plot("cluster_cooldown.txt")
#plot("cluster_heatup.txt")
#plot("heatup_zoom_cluster.txt")
#plot("cooldown_zoom_cluster.txt")
plot("00.txt")
plot("01.txt")
plot("05.txt")

#plot("cooldown.txt")

plt.show()
