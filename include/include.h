#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

// For debugging purposes
#define DEBUG 0
#define DEBUG_DUAL 0
#define DEBUG_CLAR 0

//------ DO NOT CHANGE BELOW ------
//-------- HERE BE DRAGONS --------

// number of out files
constexpr int NFILE = 2;

// information on each vertex
struct vertex {
  // vertices adjacent to it
  int adj_v[3];
  // faces it lies on
  int faces[3];
  // edges it is an endpoint of
  int edges[3];
};

// information on each face
struct face {
  // length of face (will be 5 or 6)
  int size;
  // faces it shares an edge with
  int adj_f[6];
  // vertices on face, note they will be recorded in counter clockwise order
  int vertices[6];
};

// information on each edge
struct edge {
  // vertices that form its endpoints
  // vertices[0] < vertices[1]
  int vertices[2];
};

// information on each fullerene isomer
class Fullerene {
public:
  // resize fullerene
  void Resize(int num_vertices) {
    n = num_vertices;
    // Fullerenes are 3-regular graphs, therefore a fullerene on n vertices has
    // 3n/2 edges
    num_edges = (3 * num_vertices / 2);
    // By Euler's formula, the number of faces in a fullerene is 3n/2 - n + 2
    dual_n = (num_vertices / 2 + 2);
    // resize vectors appropriately
    primal.resize(num_vertices);
    dual.resize(dual_n);
    edges.resize(num_edges);
  }
  // Attributes
  int n, dual_n, num_edges, id; // primal, dual, # of edges, and id
  vector<vertex> primal;        // planar graph information
  vector<face> dual;            // planar dual graph information
  vector<edge> edges;           // edge information
};

// information on p-anionic Clar resonance structure
class Clar_struct {
public:
  // reset values
  void Reset_vals() {
    num_res_h = 0, num_res_p = 0, num_match_e = 0;
    fill(covered_v.begin(), covered_v.end(), false);
    fill(assigned_f.begin(), assigned_f.end(), false);
  }
  // fix desired number of resonant hexagons and pentagons
  void Fix_num_res_faces(int h, int p) { res_f.resize(h + p); }
  // fix number of vertices in the graph
  void Fix_num_vert(int num_vertices) {
    // resize for number of matching edges
    match_e.resize((num_vertices - 6 * num_res_h - 5 * num_res_p) / 2);
    assigned_f.resize(num_vertices / 2 + 2); // # of faces
    covered_v.resize(num_vertices);          // # of vertices
  }
  // Attributes
  // num resonant hexagons, resonant pentagons, matching edges
  int num_res_h, num_res_p, num_match_e;
  vector<int> res_f;      // resonant faces
  vector<int> assigned_f; // # of times face has been assigned
  vector<edge> match_e;   // matching edges
  vector<int> covered_v;  // vertices covered by structure
};

// From read_and_print.cpp
void throw_error(const int n, const int h, const int p, const int graph_id,
                 string error_message);
bool read_fullerene(Fullerene(&F), const int h, const int p);
void print_primal(const int n, const vector<vertex> primal);
void print_dual(const int dual_n, const vector<face> dual);
void print_failed_test(const int graph_id, const vector<int> vec,
                       ofstream out_files_ptr[NFILE]);
void print_vec(const vector<int> vec, const string f_type);
void print_faces(const vector<int> faces, const int num_f, const string f_type);
void get_out_name(const int h, const int p, string &fname);
void open_out_file(const int h, const int p, string (&out_file_names)[NFILE],
                   ofstream out_files_ptr[NFILE]);
void close_files(ofstream out_files_ptr[NFILE]);

// From dual.cpp
int find_position(const int v, const int u, int v_adj[3]);
int counter_clockwise_walk(const int face_id, int u, int v, const int n,
                           vector<vertex>(&primal), face(&cur_face));
void construct_planar_dual(Fullerene(&F), const int h, const int p);

// From clar.cpp
void change_match(const int v_id, const int u_id, const bool match,
                  Clar_struct(&S));
bool assign_match_edges(int v_id, const Fullerene(&F), Clar_struct(&S));
bool face_term_cond_met(int *f_id, int f_count, const Fullerene(&F),
                        Clar_struct(&S), const int h, const int p,
                        ofstream out_files_ptr[NFILE]);
void change_res(const int f_id, const face f_info, bool res, Clar_struct(&S));
void assign_res_face(int *f_id, int f_count, const Fullerene(&F),
                     Clar_struct(&S), const int h, const int p,
                     ofstream out_files_ptr[NFILE]);
bool compare_face(const vector<face>(&dual), const int a, const int b);
void anionic_clar_struct_handler(const Fullerene(&F), Clar_struct(&S),
                                 const int h, const int p,
                                 ofstream out_files_ptr[NFILE]);
