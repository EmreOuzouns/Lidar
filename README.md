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

<img width="569" height="564" alt="Ekran görüntüsü 2025-12-03 161017" src="https://github.com/user-attachments/assets/e368b8ed-f791-406c-ad01-f58979e913cd" />
