#include "include.h"

int main(int argc, char *argv[]) {
  int h = atoi(argv[1]); // desired number of resonat hexagons
  int p = atoi(argv[2]); // desired number of resonat pentagons

  Fullerene F;               // isomer
  Clar_struct S;             // anionic Clar structure
  S.Fix_num_res_faces(h, p); // fixed number of resonant hex and pent

  string out_file_names[NFILE] = {"output/res_faces_", // define out files
                                  "output/graph_num_"};
  ofstream out_files_ptr[NFILE];
  open_out_file(h, p, out_file_names, out_files_ptr); // open out files

  F.id = 1;
  while (read_fullerene(F, h, p)) { // while there are isomers to read in
    construct_planar_dual(F, h, p); // construct planar dual graph
#if DEBUG
    cout << "Graph number " << graph_num << endl;
    print_primal(F.n, F.primal);
    print_dual(F.dual_n, F.dual);
#endif
    anionic_clar_struct_handler(F, S, h, p, out_files_ptr);
    F.id++;
  }
  close_files(out_files_ptr);
}
