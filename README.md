# LIDAR Data Processing Project

This project is a C implementation of a 2D LIDAR scan processing pipeline designed to turn raw polar distance measurements into structured geometric features such as clusters, wall line segments, and corner intersections. It was developed to explore how autonomous systems and robotics software filter noisy sensor scans, group physical objects, and extract reliable environmental boundaries using classical computational geometry algorithms.

The program reads raw LIDAR range data, filters out invalid or out-of-bounds readings, and converts each measurement from polar coordinates (angle and distance) into 2D Cartesian `(x, y)` points. It then groups nearby points into distinct clusters with DBSCAN, fits straight wall segments to each cluster using RANSAC, and computes valid corner intersections between extracted lines before exporting the results to output files for visualization.

---

## Technical Overview

- **Coordinate Filtering & Conversion (`main.c`):** Reads range scans, discards invalid points, and projects polar measurements into 2D Cartesian space.
- **Clustering (DBSCAN):** Separates scan points into coherent object/wall clusters while ignoring isolated noise points.
- **Line Extraction (RANSAC):** Runs `1000` RANSAC iterations per cluster (requiring a minimum of `8` inliers within `DISTANCE_THRESHOLD = 0.1`) to extract robust line segments `(A, B, C)` and endpoints.
- **Intersection & Nearest Corner Detection:** Computes intersections between extracted lines that meet at an angle of at least `60.0` degrees (`MIN_INTERSECTION_ANGLE`) and identifies the closest corner to the sensor origin.
- **Build & Run:**
  ```bash
  git clone https://github.com/EmreOuzouns/Lidar.git
  cd Lidar
  gcc main.c -o lidar -lm
  ./lidar
  ```
- **Documentation:** See `LidarProje.pdf` for the full project report.

---

## Output Visualization

The visualization below illustrates the clustered LIDAR scan points, RANSAC-fitted line segments, and detected corner intersections produced by the pipeline.

<img width="569" height="546" alt="Lidar Output" src="https://github.com/user-attachments/assets/bbd2cab0-30e1-4779-8d53-09f17623ab1d" />

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
