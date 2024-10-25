import numpy as np
import random
import itertools
import math
import matplotlib.pyplot as plt
from scipy.spatial import distance
import multiprocessing
import csv


def dist(args):
    """
    #   This function calculates the distances between randomly distributed particles on a surface.
    """
    # get the chosen random coordinated
    i = args[3]
    
    #get the first pair of random x and y coordinates    
    s1 = np.array([[args[0][i], args[1][i]]])

    #get pair of coordinates for all other random points
    s2 = np.array([[args[0][j], args[1][j]] for j in args[2] if j != i] ) 
    
    return distance.cdist(s1,s2) - d


if __name__ == '__main__':

    X = [500, 1000, 2000, 5000, 10000]#dimensions of the surface
    global_density =[0.1, 0.25, 0.5, 0.75]#% of RGD particles--

    #D = [20, 45, 90, 120]
    D = [70] #diameter of nanoparticles
    #cutoff = [50,60,70,80] #nm
    cutoff = 70 # this is what we consider nearest neighbor
    
    ratio1 = [] #Ratio 1: the RGD particles on the surface that have at least 
    # one neighboring RGD particle which is 60 nm away over THE TOTAL NUMBER OF 
    # RGD PARTICLES.
    ratio1_sd = []
    
    local_density = [1]#, 2, 4, 8] # how many RGD particle attach to 1 nanoparticle
    
    coating_average_per_local_density = [] #for each local_density, this collects the nearest neighbours 
    # per RGD particle that has atleast 1 neighbour

    for dens in local_density:
        #x, y = input('What are the dimensions of the surface (x and y, respectively)? Write in nm (tip: the bigger the surface the slower I will be, generally try 1000nm by 1000nm at max) ').split()
        x = X[2]#int(x)
        y = X[2]#int(y)
        d = D[0]#int(input('What is the diameter of your nanoparticles in nm?'))
        ncoating = 100#int(input('How many random coatings do you want to test out? (tip: the bigger the n the slower I will be!)'))
                
        coating_average = [] #how many neighbours the RGD particles have in average 
        # provided they have atleast one nearest neighbour (for n coatings)
        coating_sd = []
        
        #assuming circles side by side, get coordinates
        y_coords = []
        x_coords = []
        
        for i in range(int(d/2),int(x-d/2),d):
          for j in range(int(d/2),int(y-d/2),d):
            y_coords.append(j)
            x_coords.append(i)

        for percent in global_density:
            average_neighbour_per_rgd = [] #how many neighbours the RGD particles have in average 
            # provided they have atleast one nearest neighbour (for this coating)
            
            per_coating_ratio1 = []
            
            #loop over the number of iterations you want
            for k in range(ncoating):
                
                nn = 0 #number of RGD particles that have atleast 1 neighbour
                neighbour = 0
                #pick the random x-coordinates that are RGD particles
                random_list = random.sample(range(0,len(x_coords)), round(percent*len(x_coords)))
                #print(random_list)

                for i in range(len(random_list)):
                    
                    args = [x_coords, y_coords, random_list, random_list[i]]
                    
                    distance_ = dist(args)#calculation of neighbour local_density
                    
                    
                    for i in distance_[0]:
                        if i < 70:
                            neighbour+=1
                       
                            
                    # for higher local_density, we have to update the number of neighbours- ask Zeynep        
                    neighbour_per_local_density= dens*neighbour+(dens-1) 
                    
                    #add other particles that have RGD in nearest neighbours
                    #distance_[0] because x only has one coordinates and cdist returns the\
                    # distance in ij array form
                    if any(distance_[0] < cutoff):
                        nn += 1

                total_RGD = len(random_list)
                
                average_neighbour_per_rgd.append(neighbour_per_local_density/nn) 
                
                per_coating_ratio1.append(nn/total_RGD)
                
            '''
            Ratio 1: the RGD particles on the surface that have at least one 
            neighboring RGD particle which is 60 nm away over THE TOTAL NUMBER 
            OF RGD PARTICLES.
            '''
            ratio1.append(np.mean(per_coating_ratio1))
            
            ratio1_sd.append(np.std(per_coating_ratio1))
                
            print(str(percent) +"over!!!!")
            
            #average nn per RGD for n iterations
            coating_average.append(sum(average_neighbour_per_rgd)/ncoating)
            coating_sd.append(np.std(average_neighbour_per_rgd))
            
        coating_average_per_local_density.append(coating_average)

    fig, ax = plt.subplots(figsize =(14, 7))

    # Creating axes instance
    #ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlabel('global_density')
    ax.set_ylabel('Average number of nearest neighbours')
    # ax.set_title('Ratio of RGD particle with nearest neighbour within cutoff and total number of RGD particles ')

    # # Creating plot
    for i, value in enumerate(local_density):
        ax.errorbar(global_density,coating_average,coating_sd, label = value, marker = "o")
    ax.legend()
    plt.grid()
    plt.show()
    
    # plot ratio1
    fig, ax = plt.subplots(figsize =(7, 7))
    plt.subplots_adjust(top=0.88,
                        bottom=0.11,
                        left=0.155,
                        right=0.9,
                        hspace=0.2,
                        wspace=0.2)

    # Creating axes instance
    #ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlabel('global density', fontsize = 14)
    ax.set_ylabel(r'$\frac{RGD\,particle\,with\, nearest\, neighbour\, within\, cutoff}{Total\, RGD\, particles}$', fontsize =16)
    #ax.set_title('Ratio of RGD particle with nearest neighbour within cutoff and total number of RGD particles ')
    
    
    
    #beautifying the plot
    plt.axvline(0.1, color ='black', alpha = 0.1)
    plt.axvline(0.25, color ='black', alpha = 0.1)
    plt.axvline(0.5, color ='black', alpha = 0.1)
    # plt.axvline(0.75, color ='black', alpha = 0.1)
    
    plt.axhline(ratio1[0], color ='black', alpha = 0.1)
    plt.axhline(ratio1[1], color ='black', alpha = 0.1)
    plt.axhline(ratio1[2], color ='black', alpha = 0.1)
    # plt.axhline(ratio1[3], color ='black', alpha = 0.1)
    
    
    # Creating plot
    #for i, value in enumerate(global_density):
    ax.errorbar(global_density, ratio1, ratio1_sd, color = '#E84A5F', ls = 'dotted', marker = "o", markersize = 10)
    plt.xticks(fontsize= 14)
    plt.yticks(fontsize= 14)
    
    plt.savefig("ratio1.png", dpi = 600)
    plt.show()
    
    

    # open the file in the write mode
# with open('RGDnn_ratio_allRGD.csv', 'w') as f:
#     # create the csv writer
#     writer = csv.writer(f)

#     # write a row to the csv file
#     writer.writerow(average_nn_cutoff)
