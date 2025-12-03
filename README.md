# LIDAR Data Processing Project

This repository contains a C implementation for processing 2D LIDAR data:
- Read and filter ranges
- Convert polar to Cartesian coordinates
- Cluster using DBSCAN
- Extract lines with RANSAC
- Find intersections and nearest corner

Files:
- main.c         : C source code
- LidarProje.pdf : Project report
- README.md      : This file

## Build
gcc main.c -o lidar -lm

## Run
./lidar


<img width="569" height="546" alt="Lidar" src="https://github.com/user-attachments/assets/bbd2cab0-30e1-4779-8d53-09f17623ab1d" />
