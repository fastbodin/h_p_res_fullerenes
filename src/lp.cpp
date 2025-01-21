#include "include.h"

// int check_if_sol_valid(const Fullerene(&F), const int h, const int p) {
//   int num_res_faces = 0, res_pents = 0;
//   // for each vertex in the graph
//   for (int i = 0; i < F.n; i++) {
//     // they should be covered by the p-anionic Clar structure exactly once
//     int covered = 0;
//     for (int j = 0; j < 3; j++) {
//       // note the tolerance given to the variable assignment, i.e. > 0.99
//       // covered by matching edge
//       //if (evars[F.primal[i].edges[j]].get(GRB_DoubleAttr_X) > 0.99)
//       //  covered++;
//       // covered by resonant face
//       //if (fvars[F.primal[i].faces[j]].get(GRB_DoubleAttr_X) > 0.99)
//       //  covered++;
//     }
//     if (covered != 1) {
//       const string msg = "\nVertex " + to_string(i) + " is covered " +
//                          to_string(covered) + " times by structure.";
//       throw_error(F.n, p, F.id, msg);
//     }
//   }
//   // for each face in graph
//   for (int i = 0; i < F.dual_n; i++) {
//     //if (fvars[i].get(GRB_DoubleAttr_X) > 0.99) {
//     //  num_res_faces++;
//     //  if (F.dual[i].size == 5)
//     //    res_pents += 1;
//     //}
//   }
//   if (res_pents != p) {
//     const string msg = "\nIncorrect # of res. pents: " +
//     to_string(res_pents); throw_error(F.n, p, F.id, msg);
//   }
//   return num_res_faces;
// }
//
