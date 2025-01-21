#include "include.h"

bool pent_term(int *f_index, const Fullerene(&F), Clar_struct(&S), const int h,
               const int p, ofstream out_files_ptr[NFILE]) {
  if (p == S.num_res_p) { // correct # of resonant pentagons are assigned
    // begin assigning resonant hexagons
    return true;
  }
  if (*f_index >= 12)
    return true; // checked all pentagons
  return false;  // do not terminate
}

void assign_res_pent(int *f_index, const Fullerene(&F), Clar_struct(&S),
                     const int h, const int p, ofstream out_files_ptr[NFILE]) {
  // if the termination condition has been met for the faces
  if (pent_term(f_index, F, S, h, p, out_files_ptr))
    return;
  // check if *f_index has already been assigned (and is therefore non-res)
  if (S.p_assigned[*f_index] != 0) {
#if DEBUGCLAR
    printf("Face %d has already been assigned as non-res\n", F.pents[*f_index]);
#endif
    // move onto next pentagon
    assign_res_pent(f_index + 1, F, S, h, p, out_files_ptr);
    return;
  }
  //    // *f_index can therefore be res. We consider each case
  //    // 1) face is included in resonant set
  //    // make face resonant
  //    change_res_of_face(*f_index, F.dual[*f_index], cur_sol, F.primal, true);
  //    // move onto next pentagon
  //    assign_res_pent(f_index + 1, f_count + 1, F ,cur_sol, best_sol, hex_ids,
  //    edge_ids,
  //                    p, pents_done);
  //    // make face non-resonant
  //    change_res_of_face(*f_index, F.dual[*f_index], cur_sol, F.primal,
  //    false);
  //
  //    // 2) face is NOT included in resonant set
  //    // move onto next pentagon
  //    assign_res_pent(f_index + 1, f_count + 1, F ,cur_sol, best_sol, hex_ids,
  //    edge_ids,
  //                    p, pents_done);
  return;
}

void anionic_clar_struct_handler(const Fullerene(&F), Clar_struct(&S),
                                 const int h, const int p,
                                 ofstream out_files_ptr[NFILE]) {
  S.Fix_num_vert(F.n); // resize based on size of graph
  S.Reset_vals();      // reset values to all zero
}

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
    anionic_clar_struct_handler(F, S, h, p, out_files_ptr);

#if DEBUG
    cout << "Graph number " << graph_num << endl;
    print_primal(F.n, F.primal);
    print_dual(F.dual_n, F.dual);
    print_faces(F.pents, 12, "Pentagons: ");
    print_faces(F.hexs, F.dual_n - 12, "Hexagons: ");
#endif
    graph_num++;
  }
  close_files(out_files_ptr);
}
