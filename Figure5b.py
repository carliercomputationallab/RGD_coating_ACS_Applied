# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:37:47 2022

@author: P70073624
"""

import numpy as np
import random
import itertools
import math
import matplotlib.pyplot as plt
from scipy.spatial import distance
import multiprocessing
import pandas as pd
import csv
import os

#os.chdir(r'D:\zeynep')


def dist(args):
    """
    #   This function calculates the distances between randomly distributed particles on a surface.
    """
    # get the chosen random coordinated
    i = args[3]
    
    #get the first pair of random x and y coordinates    
    s1 = np.array([[args[0][i], args[1][i]]])

    #get pair of coordinates for all other random points
    s2 = np.array([[args[0][j], args[1][j]] for j in args[2] if j != i ]) #and abs(args[0][j] - args[0][i] ) <= cutoff and abs(args[1][j] - args[1][i] ) <= cutoff] ) 
    
    return distance.cdist(s1,s2) - d


if __name__ == '__main__':

    X = 5000#dimensions of the surface
    
    global_density =[0.25, 0.5, 1.0]#, 0.5, 1.0]#% of RGD particles--

    #D = [20, 45, 90, 120]
    d = 70 #diameter of nanoparticles
    #cutoff = [50,60,70,80] #nm
    #cutoff = 750 # 300 micro meters, this is what we consider the longest distance a cell can reach
    
    #ratio1 = [] #Ratio 1: the RGD particles on the surface that have at least 
    # one neighboring RGD particle which is 60 nm away over THE TOTAL NUMBER OF 
    # RGD PARTICLES.

    
    local_density = [1,2] # how many RGD particle attach to 1 nanoparticle
    
    eff_distance_per_global_density = []
    
    glob_05_loc_2 = []
    glob_1_loc_1 = []
    glob_025_loc_2 = []
    glob_05_loc_1 = []
    for percent in global_density:
        
        #x, y = input('What are the dimensions of the surface (x and y, respectively)? Write in nm (tip: the bigger the surface the slower I will be, generally try 1000nm by 1000nm at max) ').split()
        x = X#int(x)
        y = X#int(y)
        
        eff_distance = []
        
        ncoating = 100#int(input('How many random coatings do you want to test out? (tip: the bigger the n the slower I will be!)'))
                
        coating_average = [] #how man y neighbours the RGD particles have in average 
        # provided they have atleast one nearest neighbour (for n coatings)
        
        #assuming circles side by side, get coordinates
        y_coords = []
        x_coords = []
        
        for i in np.arange(d/2,x,d):
          for j in np.arange(d/2,y,d):
            y_coords.append(j)
            x_coords.append(i)

        for dens in local_density:
            print(dens)
            if percent == 1 and dens == 2:
               break
            avg_distance_per_coating = []
            
            #loop over the number of iterations you want
            for k in range(ncoating):
                print(k)
                avg_distance = []
                #pick the random x-coordinates that are RGD particles
                random_list = random.sample(range(0,len(x_coords)), round(percent*len(x_coords)/dens))
                #print(random_list)

                for i in range(len(random_list)):
                    
                    args = [x_coords, y_coords, random_list, random_list[i]]
                    
                    distance_ = dist(args)[0]#calculation of neighbour local_density
                    
                    avg_distance.append(dens*np.sum(distance_)/(len(random_list)))
                avg_distance_per_coating.append(np.mean(avg_distance)) # per coating!!!!   
                if percent == 0.5 and dens == 2:
                    glob_05_loc_2.append(np.mean(avg_distance))         
                elif percent == 1.0 and dens == 1:
                    glob_1_loc_1.append(np.mean(avg_distance))
                elif percent == 0.25 and dens == 2:
                    glob_025_loc_2.append(np.mean(avg_distance))
                elif percent == 0.5 and dens == 1: 
                    glob_05_loc_1.append(np.mean(avg_distance))
                    
            eff_distance.append(np.mean(avg_distance_per_coating)) # for 100 coatings average
            
                

        eff_distance_per_global_density.append(eff_distance)
        df_4_conditions = pd.DataFrame(np.transpose([glob_05_loc_2, glob_1_loc_1, glob_025_loc_2, glob_05_loc_1]))
    df_4_conditions.to_csv('df_10000.csv')    
    plt.figure()
    df_4_conditions.boxplot()
    plt.savefig('Fig5b.png')