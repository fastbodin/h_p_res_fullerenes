#include "include.h"

// For error reporting
void throw_error(const int n, const int h, const int p, const int graph_id,
                 string error_message) {
  throw runtime_error("\nError: n = " + to_string(n) + ", h = " + to_string(h) +
                      ", p = " + to_string(p) +
                      ", graph num = " + to_string(graph_id) + error_message);
}

// read in fullerene graph and populate some default values into the
// data structure
bool read_fullerene(Fullerene(&F), const int h, const int p) {
  int degree, n;
  // read in number of vertices of isomer
  if (!(cin >> n))
    return false;
  // check that number of vertices is valid
  // n should be an even number at least 20 and not equal to 22
  if (n < 20 || n == 22 || n % 2 != 0) {
    throw_error(F.n, h, p, F.id, "\nInvalid fullerene size " + to_string(n));
  }
  // resize the fullerene
  F.Resize(n);
  // for each vertex in the graph
  for (int i = 0; i < n; i++) {
    // read in degree of vertex i
    if (!(cin >> degree)) {
      throw_error(F.n, h, p, F.id,
                  "\nFailed reading vertex " + to_string(i) + "'s degree");
    }
    // check vertex degree
    if (degree != 3) {
      const string msg = "\nVertex " + to_string(i) +
                         " has invalid degree: " + to_string(degree);
      throw_error(F.n, h, p, F.id, msg);
    }
    // for each neighbour of i
    for (int j = 0; j < 3; j++) {
      // update adjacency list of vertex i
      if (!(cin >> F.primal[i].adj_v[j])) {
        const string msg = "\nFailed reading neighbour " + to_string(j) +
                           " of vertex " + to_string(i);
        throw_error(F.n, h, p, F.id, msg);
      }
      // each vertex lies on 3 faces, at this time, we do not know what their
      // ids are, we will therefore assign then as -1 to represent 'unassigned'
      F.primal[i].faces[j] = -1;
    }
  }
  return true;
}

void print_primal(const int n, const vector<vertex> primal) {
  cout << "Primal graph" << endl;
  cout << "Vert:  Neighbours         Faces           Edges" << endl;
  for (int i = 0; i < n; i++) {
    cout << setw(4) << i << ": ";
    for (int j = 0; j < 3; j++) {
      cout << setw(4) << primal[i].adj_v[j];
      if (j == 2)
        cout << " || ";
    }
    for (int j = 0; j < 3; j++) {
      cout << setw(4) << primal[i].faces[j];
      if (j == 2)
        cout << " || ";
    }
    for (int j = 0; j < 3; j++) {
      cout << setw(4) << primal[i].edges[j];
    }
    cout << endl;
  }
}

void print_dual(const int dual_n, const vector<face> dual) {
  cout << "Dual graph" << endl << "Face:       Neighbours" << endl;
  for (int i = 0; i < dual_n; i++) {
    cout << setw(4) << i << ": ";
    for (int j = 0; j < dual[i].size; j++) {
      cout << setw(4) << dual[i].adj_f[j];
    }
    cout << endl;
  }
}

void print_failed_test(const int graph_id, const vector<int> vec,
                       ofstream out_files_ptr[NFILE]) {
  out_files_ptr[1] << graph_id << endl;
  out_files_ptr[0] << vec.size();
  for (int i = 0; i < vec.size(); i++) {
    out_files_ptr[0] << " " << vec[i];
  }
  out_files_ptr[0] << endl;
}

void print_vec(const vector<int> vec, const string f_type) {
  cout << f_type;
  for (int i = 0; i < vec.size(); i++) {
    cout << " " << vec[i];
  }
  cout << endl;
}

void get_out_name(const int h, const int p, string &fname) {
  // determine first digit of h
  int h1 = h / 10;
  // determine last digit of h
  int h2 = h % 10;
  // determine first digit of p
  int p1 = p / 10;
  // determine last digit of p
  int p2 = p % 10;
  // update the file name characters
  fname[7] = (char)h1 + '0';
  fname[8] = (char)h2 + '0';
  fname[10] = (char)p1 + '0';
  fname[11] = (char)p2 + '0';
}

void open_out_file(const int h, const int p, string (&out_file_names)[NFILE],
                   ofstream out_files_ptr[NFILE]) {
  for (int i = 0; i < NFILE; i++) {
    get_out_name(h, p, out_file_names[i]);
    out_files_ptr[i].open(out_file_names[i], ios::app);
    if (!out_files_ptr[i].is_open()) {
      throw runtime_error("\nError: Could not open file " + out_file_names[i]);
    }
  }
}

void close_files(ofstream out_files_ptr[NFILE]) {
  for (int i = 0; i < NFILE; i++) {
    out_files_ptr[i].close();
  }
}
