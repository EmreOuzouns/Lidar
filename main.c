/*
   LIDAR Data Processing Project
   - read ranges
   - filtering
   - DBSCAN clustering
   - RANSAC line extraction
   - intersection point finding

   Prepared by: Emre
   GitHub-ready version.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define _USE_MATH_DEFINES 
#define MAX_RANGE 1000
#define RANSAC_ITERATIONS 1000
#define MIN_INLIERS 8 
#define DISTANCE_THRESHOLD 0.1 
#define MIN_INTERSECTION_ANGLE 60.0 

#ifdef _WIN32
    #define POPEN _popen
    #define PCLOSE _pclose
#else
    #define POPEN popen
    #define PCLOSE pclose
#endif

typedef struct {
    double m;
    double b;
    double A, B, C; 
    double start_x, start_y, end_x, end_y;
    double length;
    int cluster_id; 
    int inlier_indices[MAX_RANGE]; 
    int inlier_count;
    int original_cluster_size; 
} LineModel;

typedef struct {
    double x;
    double y;
} IntersectionPoint;

typedef struct {
    double best_distance;
    int best_i;
    int best_j;
    double best_angle;
} IntersectionResult;

void write_lines_to_file();
void write_line_points_to_file();
void write_intersections_to_file();
void save_results();

float angle_max = 0;
float angle_min = 0;
float range_max = 0;
float range_min = 0;
float angle_increment = 0;
float ranges[MAX_RANGE];
float filtered[MAX_RANGE];
float filteredAngles[MAX_RANGE];
double x_pos[MAX_RANGE];
double y_pos[MAX_RANGE];
int ranges_count = 0;
int filtered_count = 0;
int cluster_id[MAX_RANGE];
int visited[MAX_RANGE];
int found_intersection_count = 0;
int found_line_count = 0; 

LineModel found_lines[50]; 
IntersectionPoint found_intersections[100]; 

// Opens the incoming file and reads it line by line and stores to global variables

void read_parameters(FILE* fp, float* angle_min_ptr, float* angle_max_ptr, float* range_min_ptr, float* range_max_ptr, float* angle_increment_ptr) {
    char line[1000]; 
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "angle_max", 9) == 0) sscanf(line, "angle_max = %f", angle_max_ptr);

        else if (strncmp(line, "angle_min", 9) == 0) sscanf(line, "angle_min = %f", angle_min_ptr);

        else if (strncmp(line, "range_max", 9) == 0) sscanf(line, "range_max = %f", range_max_ptr);

        else if (strncmp(line, "range_min", 9) == 0) sscanf(line, "range_min = %f", range_min_ptr);

        else if (strncmp(line, "angle_increment", 15) == 0) sscanf(line, "angle_increment = %f", angle_increment_ptr);
    }
}

// When it finds the ranges line in the example file it starts reading it directly
// Uses strtok to split commas and strtof to convert to float

void read_ranges(FILE *fp, float *ranges, int *count_ptr) {
    char line[2048]; 
    int b = 0; 
    while (fgets(line, sizeof(line), fp)) {
        char *ptr = line;
        
        if (!b) {
            if (strstr(line, "ranges")) {
                b = 1; 
                ptr = strchr(line, '[');
                if (ptr) {
                    ptr++; 
                } else {
                    continue; 
                }
            } else {
                continue; 
            }
        }
        
        char *tok = strtok(ptr, ", \n\t[]");
        while (tok != NULL) {
            char *endptr;
            float value = strtof(tok, &endptr);
            
            if (endptr != tok && *tok != '\0') { 
                ranges[*count_ptr] = value;
                (*count_ptr)++;
            }
            tok = strtok(NULL, ", \n\t[]");
        }

        if (strchr(line, ']')) break;
    }
}

// Checks each point of the ranges array and if it is within min/max saves to filtered

void filter_ranges(){
    for(int i=0; i < ranges_count; i++){
        if(ranges[i] >= range_min && ranges[i] <= range_max && ranges[i] != -1.0 && !isnan(ranges[i]))
        {
            filtered[filtered_count] = ranges[i];
            filteredAngles[filtered_count] = angle_min + i * angle_increment;
            filtered_count++;
        }
    }
}

// Converts filtered distances and angles to (x,y) coordinates
// x = distance * cos(angle) and y = distance * sin(angle)

void compute_coordinates(){
    for(int a = 0; a < filtered_count; a++)
    {
       x_pos[a] = filtered[a] * cos(filteredAngles[a]);
       y_pos[a] = filtered[a] * sin(filteredAngles[a]);
    }
}

// Calculates Euclidean distance between two points from their indices

double calculate_distance(int idx1, int idx2) {
    return sqrt(pow(x_pos[idx1] - x_pos[idx2], 2) + pow(y_pos[idx1] - y_pos[idx2], 2));
}

// Finds neighbors within eps radius (region query). Returns neighbor count and fills neighbors array.

int region_query(int point_index, double eps, int neighbors[]) {
    int count = 0;
    for (int i = 0; i < filtered_count; i++) {
        if (i != point_index && calculate_distance(point_index, i) < eps) {
            neighbors[count++] = i;
        }
    }
    return count;
}

// Expands cluster for DBSCAN: if a point has enough neighbors, adds them and grows cluster

void expand_cluster(int point_index, int cid, double eps, int min_points) {
    int neighbors[MAX_RANGE]; 
    int num_neighbors =  region_query(point_index, eps, neighbors);
    int processing_queue[MAX_RANGE]; 
    int queue_count = 0;

    for(int i=0; i<num_neighbors; i++) processing_queue[queue_count++] = neighbors[i];
    int queue_idx = 0;

    while( queue_idx < queue_count) {
        int current_point_idx = processing_queue[ queue_idx++];

        if (!visited[ current_point_idx]) {
            visited[ current_point_idx] = 1;
            int new_neighbors[MAX_RANGE]; int num_new_neighbors =  region_query( current_point_idx, eps, new_neighbors);
           
            if (num_new_neighbors >= min_points) {
                for(int i = 0; i < num_new_neighbors; i++) {
                    int found = 0;
                    for(int j=0; j<queue_count; j++){ if(processing_queue[j] == new_neighbors[i]) {found = 1; break;} }
                    if(!found) processing_queue[queue_count++] = new_neighbors[i];
                }
         }
        }
        if (cluster_id[ current_point_idx] <= 0) cluster_id[ current_point_idx] = cid;
    }
}

// Iterates all points, marks noise as -1 if not enough neighbors, otherwise creates/expands cluster

int dbscan(double eps, int min_points) {
    int cid = 1;
    for (int i = 0; i < filtered_count; i++) {
        if (visited[i]) continue;
        visited[i] = 1;
        int neighbors[MAX_RANGE]; int num_neighbors =  region_query(i, eps, neighbors);
        if (num_neighbors < min_points) { cluster_id[i] = -1; }
        else {
            cluster_id[i] = cid;
            for (int j = 0; j < num_neighbors; j++) cluster_id[neighbors[j]] = cid;
            expand_cluster(i, cid, eps, min_points);
            cid++;
        }
    }
    return cid - 1;
}

// Calculates the perpendicular distance from point (Px,Py) to line Ax + By + C = 0

double point_line_distance(double Px, double Py, double A, double B, double C) {
    double num = fabs(A * Px + B * Py + C);
    double denom = sqrt(A * A + B * B);
    if (fabs(denom) < 1e-9) return 0.0;
    return num / denom;
}

// Solves analytic geometry for intersection of two lines

int find_intersection(LineModel line1, LineModel line2, IntersectionPoint* intersection) {
    double A1 = line1.A, B1 = line1.B, C1 = line1.C;
    double A2 = line2.A, B2 = line2.B, C2 = line2.C;
    double D = A1 * B2 - A2 * B1;
    if (fabs(D) < 1e-9) {
        return 0; 
    }
    intersection->x = (B1 * C2 - B2 * C1) / D;
    intersection->y = (A2 * C1 - A1 * C2) / D;
    return 1; 
}

// Finds best-fitting lines for clusters from DBSCAN using RANSAC

int ransac_find_line(int active_cluster_point_indices[], int cluster_size, int cid) {
    if (cluster_size < 2) return 0; 
    LineModel best_model = {0};
    best_model.inlier_count = 0;

    for (int i = 0; i < RANSAC_ITERATIONS; i++) {//repeat many times to find best model
        int rand_idx1 = active_cluster_point_indices[rand() % cluster_size];
        int rand_idx2 = active_cluster_point_indices[rand() % cluster_size];

        if (rand_idx1 == rand_idx2) continue;
        
        double x1 = x_pos[rand_idx1], y1 = y_pos[rand_idx1];
        double x2 = x_pos[rand_idx2], y2 = y_pos[rand_idx2];

        double A = y2 - y1;
        double B = x1 - x2;
        double C = -A * x1 - B * y1; // Ax + By + C = 0
        
        int temp_inlier_indices[MAX_RANGE];
        int temp_inlier_count = 0;

        for (int j = 0; j < cluster_size; j++) {
            int current_point_idx = active_cluster_point_indices[j];
            
            if (point_line_distance(x_pos[ current_point_idx], y_pos[ current_point_idx], A, B, C) < DISTANCE_THRESHOLD) {
                temp_inlier_indices[temp_inlier_count++] =  current_point_idx;
            }
        }
        
        if (temp_inlier_count > best_model.inlier_count) {
            best_model.A = A;
            best_model.B = B;
            best_model.C = C;

            memcpy(best_model.inlier_indices, temp_inlier_indices, temp_inlier_count * sizeof(int));
            best_model.inlier_count = temp_inlier_count;
        }
    }
    // If min inliers rule satisfied, accept the line
    if (best_model.inlier_count >= MIN_INLIERS) {
        best_model.cluster_id = cid; 
        best_model.original_cluster_size = cluster_size; 
        
        int p1_idx = -1, p2_idx = -1;
        double max_dist = -1.0;
        if (best_model.inlier_count > 1) {
            for (int i = 0; i < best_model.inlier_count; i++) {
                for (int j = i + 1; j < best_model.inlier_count; j++) {
                    double dist = calculate_distance(best_model.inlier_indices[i], best_model.inlier_indices[j]);
                    if (dist > max_dist) {
                        max_dist = dist;
                        p1_idx = best_model.inlier_indices[i];
                        p2_idx = best_model.inlier_indices[j];
             }
                }
            }
        }
        if (p1_idx != -1) {
            best_model.start_x = x_pos[p1_idx];
            best_model.start_y = y_pos[p1_idx];
            best_model.end_x = x_pos[p2_idx];
            best_model.end_y = y_pos[p2_idx];
            best_model.length = sqrt(pow(best_model.end_x - best_model.start_x, 2) + pow(best_model.end_y - best_model.start_y, 2));
        }
        
        found_lines[found_line_count++] = best_model;
        
        printf(" %d inliers found (Length: %.2fm).\n",best_model.inlier_count, best_model.length);
        
        return 1; 
    }
    
    return 0; 
}

// For plotting in the graph!

void plot_gnuplot(int total_cluster_count, double best_distance, int intersection_index_1, int intersection_index_2, double intersection_angle) {
    // Usually located here, verify:
FILE *gnuplotPipe = POPEN("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\" -persistent", "w");

    if (!gnuplotPipe) { printf("Gnuplot could not be opened!\n"); return; }
//set title : title setting
//set grid : grid setting
//set xrange: axis limit settings

    fprintf(gnuplotPipe, "set title 'LIDAR System'\n");
    fprintf(gnuplotPipe, "set xlabel 'X Position '\n");
    fprintf(gnuplotPipe, "set ylabel 'Y Position '\n");
    fprintf(gnuplotPipe, "set grid\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [-3:3]\n"); 
    fprintf(gnuplotPipe, "set yrange [-3:3]\n");
    fprintf(gnuplotPipe, "set xtics 1\n"); 
    fprintf(gnuplotPipe, "set ytics 1\n");
    fprintf(gnuplotPipe, "set key outside\n"); 

    if (found_intersection_count > 0 && best_distance > 0) {
        double label_x = found_intersections[0].x / 2.0;
        double label_y = found_intersections[0].y / 2.0;
        char label_str[50];

        sprintf(label_str, "%.2fm", best_distance);
        fprintf(gnuplotPipe, "set label 1 '%s' at %f, %f center font ',10' textcolor 'red'\n", 
                label_str, label_x, label_y);
    }else{
        fprintf(gnuplotPipe, "unset label 1\n");
    }

    if (found_intersection_count > 0 && intersection_index_1 != -1 && intersection_index_2 != -1)
     {
        char angle_label[100];
      
    } else{
        fprintf(gnuplotPipe, "unset label 2\n");
    }

    fprintf(gnuplotPipe, "plot '-' with points title 'Filtered Points' pt 7 ps 0.2 lc 'black', "
                             "'-' with points title 'Robot' pt 7 ps 2 lc 'red'");

    for (int i = 0; i < found_line_count; i++) {
        char points_title[100];
        char line_title[100];
        sprintf(points_title, "l%d points (%d)", i + 1, found_lines[i].inlier_count);
        sprintf(line_title, "l%d (%.2fm)", i + 1, found_lines[i].length);
        fprintf(gnuplotPipe, ", '-' with points title '%s' pt 7 ps 0.3 lc 'green'", points_title);
        fprintf(gnuplotPipe, ", '-' with lines title '%s' lw 2 lc 'blue'", line_title);
    }
    
    if (found_intersection_count > 0) {
        char inter_title[100];
        
        if (intersection_index_1 != -1 && intersection_index_2 != -1) {
            sprintf(inter_title, "(l%d & l%d) intersection", intersection_index_1 + 1,  intersection_index_2 + 1);  
        } else {
            sprintf(inter_title, "Nearest intersection", MIN_INTERSECTION_ANGLE);
        }

        fprintf(gnuplotPipe, ", '-' with points title '%s' pt 2 ps 2 lc 'orange'", inter_title);
        fprintf(gnuplotPipe, ", '-' with lines title 'Distance Line' lw 0.5 lc 'red' dashtype 2"); 
    }

    fprintf(gnuplotPipe, "\n");

    for (int i = 0; i < filtered_count; i++) 
    {
        fprintf(gnuplotPipe, "%f %f\n", x_pos[i], y_pos[i]);
    }

     fprintf(gnuplotPipe, "e\n");
    fprintf(gnuplotPipe, "0.0 0.0\n");
    fprintf(gnuplotPipe, "e\n");

    for (int i = 0; i < found_line_count; i++) {

        for (int j = 0; j < found_lines[i].inlier_count; j++) {
            int point_idx = found_lines[i].inlier_indices[j];
            fprintf(gnuplotPipe, "%f %f\n", x_pos[point_idx], y_pos[point_idx]);
        }
        fprintf(gnuplotPipe, "e\n");
        fprintf(gnuplotPipe, "%f %f\n%f %f\n",  
         found_lines[i].start_x, found_lines[i].start_y,
          found_lines[i].end_x, found_lines[i].end_y);
        fprintf(gnuplotPipe, "e\n"); 
    }
    
    if (found_intersection_count > 0) {
        fprintf(gnuplotPipe, "%f %f\n", found_intersections[0].x, found_intersections[0].y);
        fprintf(gnuplotPipe, "e\n");
        fprintf(gnuplotPipe, "0.0 0.0\n");

        fprintf(gnuplotPipe, "%f %f\n", found_intersections[0].x, found_intersections[0].y);
        fprintf(gnuplotPipe, "e\n");
    }

    fflush(gnuplotPipe);
    PCLOSE(gnuplotPipe);
}

int download_data(const char* url, const char* filename) {
    char command[2048];
    sprintf(command, "curl -sL \"%s\" -o \"%s\"", url, filename);
    
    printf("Downloading data: %s -> %s\n", url, filename);
    
    int ret = system(command);
    
    if (ret != 0) {
        fprintf(stderr, "Data download failed!.\n", ret);
        return 0; 
    }
    
    printf("Data successfully downloaded.\n");
    return 1; 
}

// Calls parameter reading and ranges reading to fetch data from file

int load_data(const char* filename) {
    FILE *fp_param = fopen(filename, "r");
    if (!fp_param) {
        fprintf(stderr, "'%s' file could not be opened!\n", filename);
        return 0; 
    }

    printf("Reading parameters...\n");

    read_parameters(fp_param, &angle_min, &angle_max, &range_min, &range_max, &angle_increment);
    fclose(fp_param);

    FILE *fp_range = fopen(filename, "r");
    if (!fp_range) 
    {
        fprintf(stderr, "'%s' file could not be opened!\n", filename);
        return 0; 
    }
    printf("Reading ranges...\n");

    read_ranges(fp_range, ranges, &ranges_count);
    fclose(fp_range);
    
    printf("Data reading completed. %d points read.\n", ranges_count);
    return 1; 
}

// Calls filtering and coordinate computation to make data usable
void preprocess_data() {
    printf("Filtering data and computing coordinates..\n");
    filter_ranges();
    compute_coordinates();
    printf("After filtering %d points remain.\n", filtered_count);
}

// Calls DBSCAN to split points into clusters!

int find_clusters(double epsilon, int min_points) {
    printf("Clustering in process...\n", epsilon, min_points);
    
    int total_cluster_count = dbscan(epsilon, min_points);
    
    printf("Total %d clusters found\n", total_cluster_count);
    return total_cluster_count;
}

// Analyze each cluster one by one

void extract_lines(int total_cluster_count) {
    printf("Cluster Analysis Starting\n");
    
    for (int cid = 1; cid <= total_cluster_count; cid++) {
        int active_cluster_point_indices[MAX_RANGE];
        int cluster_size = 0;
        for (int i = 0; i < filtered_count; i++) {
            if (cluster_id[i] == cid) {
                active_cluster_point_indices[cluster_size++] = i;
            }
        }
        
        printf("Cluster %d (%d points checked! ): ", cid, cluster_size);
        
        int line_found = ransac_find_line(active_cluster_point_indices, cluster_size, cid);
        
        if (!line_found) {
            printf("Not enough support found (Min: %d), not accepted as a line.\n", MIN_INLIERS);
        }
    }
    
    printf("Cluster analysis finished \n");

    int total_detected_points = 0;
    for(int i=0; i < found_line_count; i++) {
        total_detected_points += found_lines[i].inlier_count;
    }
    printf("\nTotal %d lines found.\n", found_line_count);
    printf("Total %d points assigned to these lines.\n", total_detected_points);
}

// Checks if intersection point P is on the finite segment or on infinite extension

int point_on_segment_check(IntersectionPoint P, LineModel segment, double tolerance) {
    double min_x = fmin(segment.start_x, segment.end_x) - tolerance;
    double max_x = fmax(segment.start_x, segment.end_x) + tolerance;
    double min_y = fmin(segment.start_y, segment.end_y) - tolerance;
    double max_y = fmax(segment.start_y, segment.end_y) + tolerance;

    if (P.x >= min_x && P.x <= max_x && P.y >= min_y && P.y <= max_y) {
        return 1; 
    }
    return 0; 
}

IntersectionResult intersection_analysis(){
    printf("\nIntersection Analysis Starting\n", MIN_INTERSECTION_ANGLE);
    
    found_intersection_count = 0; 
    double min_distance = 1e10; 
    int best_i = -1; 
    int best_j = -1; 
    double best_angle = 0; 

    for (int i = 0; i < found_line_count; i++) {
        for (int j = i + 1; j < found_line_count; j++) {
            IntersectionPoint temp_intersection;
           
            if (find_intersection(found_lines[i], found_lines[j], &temp_intersection)) {
                double A1 = found_lines[i].A, B1 = found_lines[i].B;
                double A2 = found_lines[j].A, B2 = found_lines[j].B;
                double dot_product = A1 * A2 + B1 * B2;
                double mag1 = sqrt(A1 * A1 + B1 * B1);
                double mag2 = sqrt(A2 * A2 + B2 * B2);
                double angle_deg = 0;

                if (mag1 > 1e-9 && mag2 > 1e-9) {
                    double cos_theta = fabs(dot_product / (mag1 * mag2));
                    if (cos_theta > 1.0) cos_theta = 1.0;
                    if (cos_theta < -1.0) cos_theta = -1.0;
                    angle_deg = acos(cos_theta) * 180.0 / M_PI;
                }// is it greater than 60 degrees!

                if (angle_deg >= MIN_INTERSECTION_ANGLE) {
                    double segment_tolerance = 5.0; 
                    
                    int on_line1 = point_on_segment_check(temp_intersection, found_lines[i], segment_tolerance);
                    int on_line2 = point_on_segment_check(temp_intersection, found_lines[j], segment_tolerance);

                    if (on_line1 && on_line2) {
                        double distance = sqrt(temp_intersection.x * temp_intersection.x + temp_intersection.y * temp_intersection.y);
                        if (distance < min_distance) {
                            min_distance = distance; 
                            best_i = i;
                            best_j = j;
                            best_angle = angle_deg; 
                            found_intersections[0] = temp_intersection; 
                            found_intersection_count = 1; 
                        }
                    } 
                   
                }
            }
        }
    }

    double gnuplot_distance_label = -1.0; 
    if (found_intersection_count > 0) {
        gnuplot_distance_label = min_distance; 
        printf("Intersection found between l%d and l%d.\n", best_i + 1, best_j + 1);
        printf("Intersection point: (%.2f, %.2f)\n", found_intersections[0].x, found_intersections[0].y);
        printf("Robot distance: %.2f m\n", min_distance);
        printf("Angle between the two lines: %.2f degrees\n", best_angle);
    } else {
        printf("No intersection point found.\n", MIN_INTERSECTION_ANGLE);
    }
    printf("Intersection Analysis Finished\n");

    IntersectionResult result = {gnuplot_distance_label, best_i, best_j, best_angle};
    return result;
}

void write_lines_to_file() {
    FILE *fp = fopen("lines.dat", "w");
    if (fp == NULL) {
        fprintf(stderr, "ERROR: lines.dat file could not be created!\n");
        return;
    }
   
    for (int i = 0; i < found_line_count; i++) {
        fprintf(fp, "%d %f %f %f %f %f %d\n",
                i + 1, 
                found_lines[i].start_x,
                found_lines[i].start_y,
                found_lines[i].end_x,
                found_lines[i].end_y,
                found_lines[i].length,
                found_lines[i].inlier_count);
    }
    fclose(fp);
}

void write_line_points_to_file() {
    FILE *fp = fopen("line_points.dat", "w");
    if (fp == NULL) {
        fprintf(stderr, "line_points.dat file could not be created!\n");
        return;
    }
    fprintf(fp, "# x y line_id\n");
    for (int i = 0; i < found_line_count; i++) {
        for (int j = 0; j < found_lines[i].inlier_count; j++) {
            int point_index = found_lines[i].inlier_indices[j];
            fprintf(fp, "%f %f %d\n",
                    x_pos[point_index],
                    y_pos[point_index],
                    i + 1); 
        }
    }
    fclose(fp);
}

void write_intersections_to_file() {
    if (found_intersection_count <= 0) {
        return;
    }
    FILE *fp = fopen("intersection.dat", "w");
    if (fp == NULL) {
        return;
    }
    fprintf(fp, "# x y\n");
    fprintf(fp, "%f %f\n",
            found_intersections[0].x,
            found_intersections[0].y);
    fclose(fp);
}

// saves results to files
void save_results() {
    write_lines_to_file();
    write_line_points_to_file();
    write_intersections_to_file();
}

void draw_results(int total_cluster_count, IntersectionResult result) {
    printf("\nDrawing graph...\n");
    plot_gnuplot(total_cluster_count, result.best_distance, result.best_i, result.best_j, result.best_angle);
}

// And the functions I call in main are here

int main() {
    srand(time(NULL));

    // fallback local file if curl fails
    const char* local_file = "scan_data_NaN.toml";

   //const char* url = " ";
   //const char* local_file = "datam.toml";

   // Using system(command) to download from url if needed
   //if (!download_data(url, local_file)) {
    //return 1; 
   //}

    // This section runs the parameter read and ranges read parts
    if (!load_data(local_file)) {
       return 1; 
    }
    
    preprocess_data();

    double epsilon = 0.07; 
    int dbscan_min_points = 2; 

    int total_clusters = find_clusters(epsilon, dbscan_min_points);

    extract_lines(total_clusters);

    IntersectionResult intersection_result = intersection_analysis();

    save_results();
    
    draw_results(total_clusters, intersection_result);

    return 0;
}
