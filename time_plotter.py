import matplotlib.pyplot as plt
import numpy as np

# nb o x o x o t0 t1....

if __name__=="__main1__":
    times = []
    file = open("results/time_mu.all","r")
    line = file.readline()
    #line.split(' ')
    while (line):
        times.append(line.split(' ')[6:])
        line = file.readline()

    
    nb_exp = len(times)
    times_mean = [0 for i in range(nb_exp)]
    for line in times:
        for i in range(len(times[0])):
            times_mean[i] += int(line[i])/nb_exp
    plt.plot(times_mean)
    plt.show()
    


if __name__=="__main__":
    times = []
    file = open("results/time_mu.all","r")
    line = file.readline()
    #line.split(' ')
    while (line):
        times.append(line.split(' ')[6:])
        line = file.readline()
        
    nb_exp = len(times)
    nb_usr = len(times[0])
    print(nb_exp,nb_usr)
    t_means = np.zeros(nb_usr)
    for i in range(nb_exp):
        t_means+=np.array(times[i],dtype='float64')
    t_means/=nb_exp
    
    plt.plot(t_means)
    plt.show()
