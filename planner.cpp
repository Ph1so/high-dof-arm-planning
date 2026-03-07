/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <algorithm>
#include <queue>
#include <map>
#include <random>

#include <tuple>
#include <string>
#include <stdexcept>
#include <regex> // For regex and split logic
#include <iostream> // cout, endl
#include <fstream> // For reading/writing files
#include <assert.h> 

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm
#define LINKLENGTH_CELLS 10

// Some potentially helpful imports
using std::vector;
using std::array;
using std::string;
using std::runtime_error;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::cout;
using std::endl;

/// @brief 
/// @param filepath 
/// @return map, x_size, y_size
tuple<double*, int, int> loadMap(string filepath) {
	std::FILE *f = fopen(filepath.c_str(), "r");
	if (f) {
	}
	else {
		printf("Opening file failed! \n");
		throw runtime_error("Opening map file failed!");
	}
	int height, width;
	if (fscanf(f, "height %d\nwidth %d\n", &height, &width) != 2) {
		throw runtime_error("Invalid loadMap parsing map metadata");
	}
	
	////// Go through file and add to m_occupancy
	double* map = new double[height*width];

	double cx, cy, cz;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			char c;
			do {
				if (fscanf(f, "%c", &c) != 1) {
					throw runtime_error("Invalid parsing individual map data");
				}
			} while (isspace(c));
			if (!(c == '0')) { 
				map[y+x*width] = 1; // Note transposed from visual
			} else {
				map[y+x*width] = 0;
			}
		}
	}
	fclose(f);
	return make_tuple(map, width, height);
}

// Splits string based on deliminator
vector<string> split(const string& str, const string& delim) {   
		// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
		const std::regex ws_re(delim);
		return { std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator() };
}


double* doubleArrayFromString(string str) {
	vector<string> vals = split(str, ",");
	double* ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i) {
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double* v1, double *v2, int size) {
    for (int i = 0; i < size; ++i) {
        if (abs(v1[i]-v2[i]) > 1e-3) {
            cout << endl;
            return false;
        }
    }
    return true;
}

typedef struct {
	int X1, Y1;
	int X2, Y2;
	int Increment;
	int UsingYIndex;
	int DeltaX, DeltaY;
	int DTerm;
	int IncrE, IncrNE;
	int XIndex, YIndex;
	int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size) {
	double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params) {
	params->UsingYIndex = 0;

	if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
		(params->UsingYIndex)++;

	if (params->UsingYIndex)
		{
			params->Y1=p1x;
			params->X1=p1y;
			params->Y2=p2x;
			params->X2=p2y;
		}
	else
		{
			params->X1=p1x;
			params->Y1=p1y;
			params->X2=p2x;
			params->Y2=p2y;
		}

	 if ((p2x - p1x) * (p2y - p1y) < 0)
		{
			params->Flipped = 1;
			params->Y1 = -params->Y1;
			params->Y2 = -params->Y2;
		}
	else
		params->Flipped = 0;

	if (params->X2 > params->X1)
		params->Increment = 1;
	else
		params->Increment = -1;

	params->DeltaX=params->X2-params->X1;
	params->DeltaY=params->Y2-params->Y1;

	params->IncrE=2*params->DeltaY*params->Increment;
	params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
	params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

	params->XIndex = params->X1;
	params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y) {
	if (params->UsingYIndex) {
        *y = params->XIndex;
        *x = params->YIndex;
        if (params->Flipped)
            *x = -*x;
    }
	else {
        *x = params->XIndex;
        *y = params->YIndex;
        if (params->Flipped)
            *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params) {
	if (params->XIndex == params->X2) {
        return 0;
    }
	params->XIndex += params->Increment;
	if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
		params->DTerm += params->IncrE;
	else {
        params->DTerm += params->IncrNE;
        params->YIndex += params->Increment;
	}
	return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
			 int x_size, int y_size) {
	bresenham_param_t params;
	int nX, nY; 
	short unsigned int nX0, nY0, nX1, nY1;

	//printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
		
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

	//printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
			return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
			 int x_size, int y_size) {
    double x0,y0,x1,y1;
    int i;
		
	 //iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
	y1 = 0;
	for(i = 0; i < numofDOFs; i++){
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
			return 0;
	}    
	return 1;
}

double calc_cost(double** plan, int numofDOFs, int planlength)
{
	double total_cost = 0;
	for (int i = 1; i < planlength; i++)
	{
		for (int n = 0; n < numofDOFs; n++)
		{
			total_cost += fabs(plan[i][n] - plan[i-1][n]);
		}
	}
	return total_cost;
}


static void linear_interp_planner(
			double* map,
			int x_size,
			int y_size,
			double* armstart_anglesV_rad,
			double* armgoal_anglesV_rad,
            int numofDOFs,
            double*** plan,
            int* planlength)
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
		
    //for now just do straight interpolation between start and goal checking for the validity of samples

    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("The arm is already at the goal\n");
        return;
    }
	int countNumInvalid = 0;
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size)) {
			++countNumInvalid;
        }
    }
	printf("Linear interpolation collided at %d instances across the path\n", countNumInvalid);
    *planlength = numofsamples;
    
    return;
}

bool linear_interp(double* map, int x_size, int y_size, double* start, double* end, int numofDOFs) {
    std::vector<double> diffs(numofDOFs);
    double max_diff = 0;

    for (int j = 0; j < numofDOFs; j++) {
        double d = end[j] - start[j];
        
        // wrap the difference to [-PI, PI] - finds the shortest way around the circle
        while (d > PI)  d -= 2 * PI;
        while (d < -PI) d += 2 * PI;
        
        diffs[j] = d;
        if (fabs(d) > max_diff) max_diff = fabs(d);
    }
    
    // Use max_diff to determine sampling density
    int numofsamples = (int)(max_diff / (PI / 20));
    if (numofsamples < 2) return true;

    std::vector<double> temp_config(numofDOFs);
    for (int i = 0; i < numofsamples; i++) {
        double ratio = (double)i / (numofsamples - 1);
        for (int j = 0; j < numofDOFs; j++) {
            // interpolate using the shortest path difference
            double val = start[j] + ratio * diffs[j];
            
            // keep the resulting angle within [0, 2*PI] for the collision checker
            while (val < 0) val += 2 * PI;
            while (val >= 2 * PI) val -= 2 * PI;
            
            temp_config[j] = val;
        }
        if (!IsValidArmConfiguration(temp_config.data(), numofDOFs, map, x_size, y_size)) {
            return false;
        }
    }
    return true;
}

/** PRM Planner
 * 
 * NOTE: in progress
*/
static void prm(
			double* map,
			int x_size,
			int y_size,
			double* armstart_anglesV_rad,
			double* armgoal_anglesV_rad,
            int numofDOFs,
            double*** plan,
            int* planlength,
			int n,
			int k
        	)
{
	// printf("Running PRM\n");
	// auto start_time = std::chrono::high_resolution_clock::now();

	// no plan by default
	*plan = NULL;
	*planlength = 0;

	// constants
	typedef double *Node;
    using State = std::pair<double, int>;  
	
	const string mapfile = "map2.txt";
	// const int n = 10000;
	// const int k = 200;
	const long max_rand = 1000000L;
	const double lower_bound = 0;
	const double upper_bound = PI;

	// init vars
	double distance;
	srandom(time(NULL));

	// some helpers
	auto is_connected = [&](Node start, Node end) -> bool {
		return linear_interp(map, x_size, y_size, start, end, numofDOFs);
    };
	
	// init graph
	std::map<Node, int> mappings; 
	Node *graph_set = (Node *)calloc(sizeof(Node), n);
	double **graph_matrix = (double**)malloc(sizeof(double*)*n); // alloc rows
	for (int i = 0; i < n; i++)
	{
		graph_matrix[i] = (double*)calloc(sizeof(double), n); // alloc columns
	}

	// is it faster to use std::vector? e.g. std::vector<Node **> gm[n];

	// setup random number engines 
	std::default_random_engine generator(time(NULL));
	std::uniform_real_distribution<double> uniform(0.0, 2 * PI);
	// Standard deviation (sigma) determines how close nodes are to obstacles
	// A sigma of PI/10 (18 degrees) is a good starting point
	std::normal_distribution<double> gaussian(0.0, PI / 10);

	// node generation
	graph_set[0] = armstart_anglesV_rad; // Start node
	graph_set[1] = armgoal_anglesV_rad;  // Goal node
	// Start your random loop from i = 2
	for (int i = 2; i < n; i++) 
	{
		Node q1 = (Node)malloc(numofDOFs * sizeof(double));
		Node q2 = (Node)malloc(numofDOFs * sizeof(double));
		
		// Generate q1 (Uniform)
		for (int j = 0; j < numofDOFs; j++) q1[j] = uniform(generator);
		
		// Generate q2 (Gaussian around q1)
		for (int j = 0; j < numofDOFs; j++) {
			q2[j] = q1[j] + gaussian(generator);
			// Wrap to [0, 2PI]
			while (q2[j] < 0) q2[j] += 2 * PI;
			while (q2[j] >= 2 * PI) q2[j] -= 2 * PI;
		}

		bool q1_valid = IsValidArmConfiguration(q1, numofDOFs, map, x_size, y_size);
		bool q2_valid = IsValidArmConfiguration(q2, numofDOFs, map, x_size, y_size);

		Node node = NULL;
		// gaussian sampling: only care if exactly one is in collision
		if (q1_valid && !q2_valid) {
			node = q1;
			free(q2); 
		} else if (!q1_valid && q2_valid) {
			node = q2;
			free(q1); 
		} else {
			// Both valid or both invalid: discard both
			free(q1);
			free(q2);
			i--; // Decrement to try again
			continue;
		}

		//add to graph
		graph_set[i] = node;
		mappings[node] = i;
	}

	// connect the nodes
	for (int i = 0; i < n; i++)
	{
		// check if node can be connected to the k closest nodes by doing lin interp 
		// first, get distances from neighbors
    	std::priority_queue<State, std::vector<State>, std::greater<State>> distances;
		Node node = graph_set[i];
		for (int neighbor_i = 0; neighbor_i < i; neighbor_i++)
		{
			// remeber to ignore itself
			// if (neighbor_i == i) continue; 

			Node neighbor = graph_set[neighbor_i];
			if (neighbor == NULL) continue;

			distance = 0;
			for (int angle_i = 0; angle_i < numofDOFs; angle_i++)
			{
				distance += (neighbor[angle_i] - node[angle_i]) * (neighbor[angle_i] - node[angle_i]);
			}
			distances.push({distance, neighbor_i});
		}

		// second, check k closest neighbors
		for (int k_i = 0; k_i < k && !distances.empty(); k_i++) {
			int k_candidate = distances.top().second;
			double dist = distances.top().first;
			distances.pop();

			if (is_connected(graph_set[k_candidate], node)) 
			{
				graph_matrix[k_candidate][i] = dist;
				graph_matrix[i][k_candidate] = dist;
			}
		}
	}

	// save PRM graph for visualization (not necessary for actual planning)
	// string map_basename = mapfile.substr(mapfile.find_last_of("/\\") + 1);
	// string prm_filename = "prm_graph_" + map_basename;
	// std::ofstream prm_file(prm_filename, std::ios::trunc);
	// if (prm_file.is_open()) {
	// 	// Count valid nodes
	// 	int valid_count = 0;
	// 	for (int i = 0; i < n; i++) {
	// 		if (graph_set[i] != nullptr) valid_count++;
	// 	}

	// 	prm_file << "DOFS " << numofDOFs << "\n";
	// 	prm_file << "NODES " << valid_count << "\n";
	// 	for (int i = 0; i < n; i++) {
	// 		if (graph_set[i] == nullptr) continue;
	// 		prm_file << i;
	// 		for (int j = 0; j < numofDOFs; j++) {
	// 			prm_file << " " << graph_set[i][j];
	// 		}
	// 		prm_file << "\n";
	// 	}

	// 	prm_file << "EDGES\n";
	// 	for (int i = 0; i < n; i++) {
	// 		if (graph_set[i] == nullptr) continue;
	// 		for (int j = i + 1; j < n; j++) {
	// 			if (graph_set[j] == nullptr) continue;
	// 			if (graph_matrix[i][j] != 0) {
	// 				prm_file << i << " " << j << "\n";
	// 			}
	// 		}
	// 	}
	// 	prm_file.close();
	// 	cout << "PRM graph saved to " << prm_filename << endl;
	// }

	// A* search
    std::vector<double> g_score(n, std::numeric_limits<double>::infinity());
    std::vector<int> parent(n, -1);
    
    // priority queue stores pair<f_score, node_index>
    using AStarState = std::pair<double, int>;
    std::priority_queue<AStarState, std::vector<AStarState>, std::greater<AStarState>> open_set;

    int start_node = 0;
    int goal_node = 1;

    g_score[start_node] = 0;
    
    auto get_heuristic = [&](int idx) {
        double h = 0;
        for (int d = 0; d < numofDOFs; d++) {
            double diff = graph_set[idx][d] - graph_set[goal_node][d];
            h += diff * diff;
        }
        return sqrt(h);
    };

    open_set.push({get_heuristic(start_node), start_node});

    bool found_path = false;
    while (!open_set.empty()) {
        int current = open_set.top().second;
        open_set.pop();

        if (current == goal_node) {
            found_path = true;
            break;
        }

        for (int neighbor = 0; neighbor < n; neighbor++) {
            double edge_weight = graph_matrix[current][neighbor];
            if (edge_weight > 0) { // If there's an edge
                double tentative_g = g_score[current] + sqrt(edge_weight);

                if (tentative_g < g_score[neighbor]) {
                    parent[neighbor] = current;
                    g_score[neighbor] = tentative_g;
                    double f_score = tentative_g + get_heuristic(neighbor);
                    open_set.push({f_score, neighbor});
                }
            }
        }
    }

    // path reconstruction with interp
    if (found_path) {
        std::vector<int> path_indices;
        int curr = goal_node;
        while (curr != -1) {
            path_indices.push_back(curr);
            curr = parent[curr];
        }
        std::reverse(path_indices.begin(), path_indices.end());

        // max 16 degree change per step
        const double MAX_STEP = 8.0 * (PI / 180.0);
        std::vector<std::vector<double>> interpolated_path;

        for (size_t i = 0; i < path_indices.size() - 1; i++) {
            double* start_node = graph_set[path_indices[i]];
            double* end_node = graph_set[path_indices[i+1]];

            // maximum angular distance between any joint
            double max_dist = 0;
            for (int d = 0; d < numofDOFs; d++) {
                double diff = fabs(end_node[d] - start_node[d]);
                if (diff > max_dist) max_dist = diff;
            }

            // how many segments to satisfy the 16 degree limit
            int steps = (int)ceil(max_dist / MAX_STEP);
            if (steps < 1) steps = 1;

            // interp (excluding the last point, which is the start of the next segment)
            for (int s = 0; s < steps; s++) {
                std::vector<double> state(numofDOFs);
                double t = (double)s / steps;
                for (int d = 0; d < numofDOFs; d++) {
                    state[d] = start_node[d] + t * (end_node[d] - start_node[d]);
                }
                interpolated_path.push_back(state);
            }
        }
        
        // add the goal node
        std::vector<double> goal_state(numofDOFs);
        for(int d = 0; d < numofDOFs; d++) goal_state[d] = graph_set[goal_node][d];
        interpolated_path.push_back(goal_state);

        // convert std::vector to the double** format
        *planlength = (int)interpolated_path.size();
        *plan = (double**)malloc((*planlength) * sizeof(double*));

        for (int i = 0; i < *planlength; i++) {
            (*plan)[i] = (double*)malloc(numofDOFs * sizeof(double));
            for (int d = 0; d < numofDOFs; d++) {
                (*plan)[i][d] = interpolated_path[i][d];
            }
        }
        // printf("A* found path. After 16-deg interpolation, plan length is %d\n", *planlength);
    } else {
        printf("A* failed to find a path.\n");
    }

	// auto end_time = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> elapsed = end_time - start_time;
    // printf("time: %.4f ms\n", elapsed.count());

    return;
}

#include <chrono>

struct OptResult {
    int n;
    int k;
    double cost;
    double time;
};

// A wrapper for your PRM that returns performance metrics
OptResult evaluate_prm(int n, int k, double* map, int x_size, int y_size, 
                      double* start, double* goal, int dofs) {
    double** plan = NULL;
    int planlength = 0;
    
    auto start_t = std::chrono::high_resolution_clock::now();
    
    // Call your prm function (modified to take n and k)
    prm(map, x_size, y_size, start, goal, dofs, &plan, &planlength, n, k);
    
    auto end_t = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_t - start_t;

    double cost = (planlength > 0) ? calc_cost(plan, dofs, planlength) : 1e9;
    
    // Clean up plan memory to prevent leaks during optimization
    if(plan) {
        for(int i=0; i<planlength; i++) free(plan[i]);
        free(plan);
    }

    return {n, k, cost, elapsed.count()};
}

void optimize_planner(double* map, int x_size, int y_size, double* start,
                      double* goal, int dofs, double*** final_plan, int* final_length) {

    // === Bayesian Optimization over (n, k) ===
    // Surrogate: Gaussian Process with RBF kernel
    // Acquisition: Expected Improvement (minimization)
    const double N_MIN = 1000, N_MAX = 8000;
    const double K_MIN = 10,   K_MAX = 200;
    const double TIME_LIMIT = 4.5; // seconds (buffer below 5000ms limit)
    const int INIT_SAMPLES = 3;
    const int MAX_ITER = 500;
    const int RUNS_PER_EVAL = 5;
    const double GP_NOISE = 1e-4;
    const double LENGTH_SCALE = 0.3; // RBF kernel length scale in normalized space

    // Normalize helpers
    auto norm_n = [&](double n) { return (n - N_MIN) / (N_MAX - N_MIN); };
    auto norm_k = [&](double k) { return (k - K_MIN) / (K_MAX - K_MIN); };

    // RBF kernel between two normalized points
    auto rbf = [&](double n1, double k1, double n2, double k2) {
        double dn = n1 - n2, dk = k1 - k2;
        return std::exp(-(dn*dn + dk*dk) / (2.0 * LENGTH_SCALE * LENGTH_SCALE));
    };

    std::default_random_engine gen(time(NULL));
    std::uniform_real_distribution<double> rand_n(N_MIN, N_MAX);
    std::uniform_real_distribution<double> rand_k(K_MIN, K_MAX);

    std::vector<double> obs_n, obs_k, obs_cost;
    OptResult best_found = {0, 0, 1e9, 0.0};

    auto add_observation = [&](int test_n, int test_k) {
        double total_cost = 0.0, total_time = 0.0;
        for (int r = 0; r < RUNS_PER_EVAL; r++) {
            OptResult res = evaluate_prm(test_n, test_k, map, x_size, y_size, start, goal, dofs);
            total_cost += res.cost;
            total_time += res.time;
        }
        double avg_cost = total_cost / RUNS_PER_EVAL;
        double avg_time = total_time / RUNS_PER_EVAL;

        obs_n.push_back(norm_n(test_n));
        obs_k.push_back(norm_k(test_k));
        // Penalize infeasible (too slow) runs with a large cost
        double penalized = (avg_time < TIME_LIMIT) ? avg_cost : 1e9;
        obs_cost.push_back(penalized);
        if (avg_time < TIME_LIMIT && avg_cost < best_found.cost) {
            best_found = {test_n, test_k, avg_cost, avg_time};
            printf("[BO] New Best: n=%d, k=%d, AvgCost=%.2f, AvgTime=%.2fs\n",
                   test_n, test_k, avg_cost, avg_time);
        }
    };

    // --- Phase 1: random initial samples ---
    for (int i = 0; i < INIT_SAMPLES; i++)
        add_observation((int)rand_n(gen), (int)rand_k(gen));

    // --- Phase 2: BO iterations ---
    for (int iter = INIT_SAMPLES; iter < MAX_ITER; iter++) {
		printf("%d\n", iter);
        int m = (int)obs_n.size();

        // Build K matrix (m x m) with noise on diagonal
        std::vector<std::vector<double>> L(m, std::vector<double>(m, 0.0));
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                L[i][j] = rbf(obs_n[i], obs_k[i], obs_n[j], obs_k[j]) + (i == j ? GP_NOISE : 0.0);

        // Cholesky: L overwritten with lower triangular factor
        for (int i = 0; i < m; i++) {
            for (int j = 0; j <= i; j++) {
                double s = L[i][j];
                for (int p = 0; p < j; p++) s -= L[i][p] * L[j][p];
                L[i][j] = (i == j) ? std::sqrt(std::max(s, 1e-12)) : s / L[j][j];
            }
            for (int j = i + 1; j < m; j++) L[i][j] = 0.0;
        }

        // Solve K * alpha = obs_cost via Cholesky (forward then back sub)
        std::vector<double> tmp(m), alpha(m);
        for (int i = 0; i < m; i++) {
            double s = obs_cost[i];
            for (int j = 0; j < i; j++) s -= L[i][j] * tmp[j];
            tmp[i] = s / L[i][i];
        }
        for (int i = m - 1; i >= 0; i--) {
            double s = tmp[i];
            for (int j = i + 1; j < m; j++) s -= L[j][i] * alpha[j];
            alpha[i] = s / L[i][i];
        }

        double f_best = *std::min_element(obs_cost.begin(), obs_cost.end());

        // Grid search for max Expected Improvement
        const int GRID = 35;
        double best_ei = -1.0;
        double cand_n = obs_n[0], cand_k = obs_k[0];

        for (int gi = 0; gi < GRID; gi++) {
            for (int gj = 0; gj < GRID; gj++) {
                double cn = gi / (double)(GRID - 1);
                double ck = gj / (double)(GRID - 1);

                // GP posterior mean
                double mu = 0.0;
                for (int i = 0; i < m; i++)
                    mu += rbf(cn, ck, obs_n[i], obs_k[i]) * alpha[i];

                // GP posterior variance via v = L^{-1} k_star
                std::vector<double> v(m);
                for (int i = 0; i < m; i++) {
                    double s = rbf(cn, ck, obs_n[i], obs_k[i]);
                    for (int j = 0; j < i; j++) s -= L[i][j] * v[j];
                    v[i] = s / L[i][i];
                }
                double var = 1.0; // k(x*, x*) = 1 for normalized RBF
                for (int i = 0; i < m; i++) var -= v[i] * v[i];
                double sigma = std::sqrt(std::max(var, 0.0));

                // Expected Improvement (minimization form)
                double ei = 0.0;
                if (sigma > 1e-10) {
                    double z = (f_best - mu) / sigma;
                    double phi = std::exp(-0.5 * z * z) / std::sqrt(2.0 * M_PI);
                    double Phi = 0.5 * (1.0 + std::erf(z / std::sqrt(2.0)));
                    ei = (f_best - mu) * Phi + sigma * phi;
                }
                if (ei > best_ei) { best_ei = ei; cand_n = cn; cand_k = ck; }
            }
        }

        // Denormalize candidate
        int test_n = std::max((int)N_MIN, std::min((int)N_MAX,
                        (int)(cand_n * (N_MAX - N_MIN) + N_MIN)));
        int test_k = std::max((int)K_MIN, std::min((int)K_MAX,
                        (int)(cand_k * (K_MAX - K_MIN) + K_MIN)));
        add_observation(test_n, test_k);
    }

    // Final run with best found parameters
    if (best_found.n > 0)
        prm(map, x_size, y_size, start, goal, dofs, final_plan, final_length,
            best_found.n, best_found.k);
    else
        prm(map, x_size, y_size, start, goal, dofs, final_plan, final_length, 5000, 20);
	printf("Best parameters n: %d, and k: %d with cost: %f\n", best_found.n, best_found.k, calc_cost(*final_plan, dofs, *final_length));
}

/** Your final solution will be graded by an grading script which will
 * send the default 6 arguments:
 *    map, numOfDOFs, commaSeparatedStartPos, commaSeparatedGoalPos, 
 *    whichPlanner, outputFilePath
 * An example run after compiling and getting the planner.out executable
 * >> ./planner.out map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 0 output.txt
 * See the hw handout for full information.
 * If you modify this for testing (e.g. to try out different hyper-parameters),
 * make sure it can run with the original 6 commands.
 * Programs that do not will automatically get a 0.
 * */
int main(int argc, char** argv) {
	double* map;
	int x_size, y_size;

	tie(map, x_size, y_size) = loadMap(argv[1]);
	const int numOfDOFs = std::stoi(argv[2]);
	double* startPos = doubleArrayFromString(argv[3]);
	double* goalPos = doubleArrayFromString(argv[4]);
	int whichPlanner = std::stoi(argv[5]);
	string outputFile = argv[6];

	if(!IsValidArmConfiguration(startPos, numOfDOFs, map, x_size, y_size)||
			!IsValidArmConfiguration(goalPos, numOfDOFs, map, x_size, y_size)) {
		throw runtime_error("Invalid start or goal configuration!\n");
	}

	///////////////////////////////////////
	//// Feel free to modify anything below. Be careful modifying anything above.

	double** plan = NULL;
	int planlength = 0;
	auto start_time = std::chrono::high_resolution_clock::now();
	if (whichPlanner == 3) {
		printf("Running PRM\n");
		prm(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength, 2000, 160);
	}
	else linear_interp_planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
	auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end_time - start_time;
    printf("time: %.4f ms\n", elapsed.count());
	printf("cost: %f\n", calc_cost(plan, numOfDOFs, planlength));

	// optimize_planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);

	//// Feel free to modify anything above.
	//// If you modify something below, please change it back afterwards as the 
	//// grading script will not work.
	///////////////////////////////////////

    // Your solution's path should start with startPos and end with goalPos
    if (!equalDoubleArrays(plan[0], startPos, numOfDOFs) || 
    	!equalDoubleArrays(plan[planlength-1], goalPos, numOfDOFs)) {
		throw std::runtime_error("Start or goal position not matching");
	}

	/** Saves the solution to output file
	 * Do not modify the output log file output format as it is required for visualization
	 * and for grading.
	 */
	std::ofstream m_log_fstream;
	m_log_fstream.open(outputFile, std::ios::trunc); // Creates new or replaces existing file
	if (!m_log_fstream.is_open()) {
		throw std::runtime_error("Cannot open file");
	}
	m_log_fstream << argv[1] << endl; // Write out map name first
	/// Then write out all the joint angles in the plan sequentially
	for (int i = 0; i < planlength; ++i) {
		for (int k = 0; k < numOfDOFs; ++k) {
			m_log_fstream << plan[i][k] << ",";
		}
		m_log_fstream << endl;
	}
}
