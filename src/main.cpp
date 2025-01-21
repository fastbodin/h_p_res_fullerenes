#include "include.h"

int main(int argc, char *argv[]) {
  int h = atoi(argv[1]); // desired number of resonat hexagons
  int p = atoi(argv[2]); // desired number of resonat pentagons

  Fullerene F;               // isomer
  Clar_struct S;             // anionic Clar structure
  S.Fix_num_res_faces(h, p); // fixed number of resonant hex and pent

  string out_file_names[NFILE] = {"output/hh_pp_res_faces", // define out files
                                  "output/hh_pp_graph_num"};
  ofstream out_files_ptr[NFILE];
  open_out_file(h, p, out_file_names, out_files_ptr); // open out files

  int graph_num = 0;
  while (read_fullerene(F, h, p)) { // while there are isomers to read in
    F.id = graph_num;
    construct_planar_dual(F, h, p); // construct planar dual graph
#if DEBUG
    cout << "Graph number " << graph_num << endl;
    print_primal(F.n, F.primal);
    print_dual(F.dual_n, F.dual);
#endif

    anionic_clar_struct_handler(F, S, h, p, out_files_ptr);
    graph_num++;
  }
  close_files(out_files_ptr);
}
