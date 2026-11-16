# LIDAR Data Processing Project

A C implementation for processing 2D LIDAR scan data, clustering points, extracting wall lines, and finding intersection corners.

<img width="569" height="546" alt="Lidar Output" src="https://github.com/user-attachments/assets/bbd2cab0-30e1-4779-8d53-09f17623ab1d" />

## What It Does

- Reads raw 2D LIDAR range measurements and filters invalid readings.
- Converts polar coordinates (angle, distance) into Cartesian `(x, y)` coordinates.
- Groups nearby points into clusters using DBSCAN.
- Fits line segments to each cluster using RANSAC (`1000` iterations, minimum `8` inliers).
- Computes line intersection points (filtering by minimum intersection angle of `60` degrees) and identifies the nearest corner.
- Exports processed line segments, cluster points, and intersection coordinates to output files for plotting.

## Files

- `main.c` - Full C implementation of the LIDAR processing pipeline
- `LidarProje.pdf` - Project report
- `README.md` - This file

## Build & Run

```bash
gcc main.c -o lidar -lm
./lidar
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
