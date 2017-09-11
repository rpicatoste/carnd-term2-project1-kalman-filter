#**Self-Driving Car Nanodegree Program - Term 2** 
##**Project 1 - Extended Kalman Filter** 
**Ricardo Picatoste**

## Notes
The project has been compiled in Windows 10 using the "Bash on Ubuntu on Windows", generating the make file with "cmake CMakeLists.txt" and then running make. 

I have used a tab width of 4 spaces, taking care of the matrices alignment. I hope it's ok, I prefer it like that for programming.


## Results

The program runs with the Dataset 1 and 2, even restarting them without closing the simulator or the ExtendedKF executable, as it will re-initialize matrices and states as needed.

The RMSE achieved is below the required maximum limit. To achieve these values, it has to run for a good part of the dataset (where there are just a few points, the RMSE is higher).

Below there is a capture after finishing dataset1:
![alt text](./results/after_dataset1.png "After Dataset1")


Below there is a capture after finishing dataset2:
![alt text](./results/after_dataset2.png "After Dataset2")