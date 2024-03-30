This folder contains the code for the paper ”Clustering Computer Mouse Tracking Data with Informed Hierarchical Shrinkage Partition Priors" by Ziyi Song, Weining Shen, Marina Vannucci, Alexandria Baldizon, Paul M. Cinciripini, Francesco Versace, and Michele Guindani (2023).

Feel free to download the zip file. Codes used for the paper should be in it.

The folders Simulation 1a, Simulation 1b, and Simulation 2 correspond to the simulation scenarios in the paper. Each folder contains the needed codes and simulated data and results. You can also see a Description text file in each folder that explains details about the folder.


The folder "Functions" contains the main functions implemented for our HSP method. Followings are brief explanations about the functions.

"SimData.R": Data simulation. It might be different across simulation scenarios. See scenarios in the paper. 
"relabel.R": Given labels of a vector of items, it relabel the items using incremental label numbers, while the clustering unchanged. For example, {3,3,1,2,2} --> {1,1,2,3,3}.


