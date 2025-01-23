#include "include.h"

int main(int argc, char *argv[]) {
  const int h = 0; // desired number of resonat hexagons
  const int p = 2; // desired number of resonat pentagons

  Fullerene F;                    // isomer
  read_fullerene(F, h, p);        // read in isomer C48:32
  construct_planar_dual(F, h, p); // construct planar dual graph

  Clar_struct S;             // anionic Clar structure
  S.Fix_num_res_faces(h, p); // fixed number of resonant hex and pent
  S.Fix_num_vert(F.n);       // resize based on size of graph
  S.Reset_vals();            // reset values

  // isomer C48:32 is not 2-anionic resonant
  // it has 3 patches of pentagons (one of which is isolated: face 11)
  // the patches are:
  int pent_patch_1[5] = {0, 1, 2, 3, 4};
  int pent_patch_2[6] = {18, 19, 21, 22, 24, 25};
  int p_1, p_2;
  // the deletion of the vertices from any pair of pentagons,
  // one from patch_1 and one from patch_2, results in a a graph
  // that contains a perfect matching
  for (int i = 0; i < 5; i++) {
    p_1 = pent_patch_1[i];
    change_res(p_1, F.dual[p_1], true, S);
    for (int j = 0; j < 6; j++) {
      p_2 = pent_patch_2[j];
      change_res(p_2, F.dual[p_2], true, S);
      if (!assign_match_edges(0, F, S)) {
        throw_error(F.n, h, p, 1, "Perfect matching not found");
      }
      change_res(p_2, F.dual[p_2], false, S);
    }
    change_res(p_1, F.dual[p_1], false, S);
  }
  // since patch_1 does not contain two independent pentagons,
  // we consider choices of two independent pentagons within patch_2
  // one of the face from patch_2 must be: 18, 19, or 24
  // chooose 18 or 19
  for (int i = 0; i < 2; i++) {
    p_1 = pent_patch_2[i];
    change_res(p_1, F.dual[p_1], true, S);
    // if 18 or 19 is chosen from patch_2, then there are
    // three choices for the second pentagon
    change_res(21, F.dual[21], true, S);
    if (!assign_match_edges(0, F, S)) {
      throw_error(F.n, h, p, 1, "Perfect matching not found");
    }
    change_res(21, F.dual[21], false, S);

    change_res(22, F.dual[22], true, S);
    if (!assign_match_edges(0, F, S)) {
      throw_error(F.n, h, p, 1, "Perfect matching not found");
    }
    change_res(22, F.dual[22], false, S);

    change_res(25, F.dual[25], true, S);
    if (!assign_match_edges(0, F, S)) {
      throw_error(F.n, h, p, 1, "Perfect matching not found");
    }
    change_res(25, F.dual[25], false, S);

    change_res(p_1, F.dual[p_1], false, S);
  }
  // choose 24
  change_res(24, F.dual[24], true, S);
  // if 24 is chosen from patch_2, then there are
  // two choices for the second pentagon
  change_res(21, F.dual[21], true, S);
  if (!assign_match_edges(0, F, S)) {
    throw_error(F.n, h, p, 1, "Perfect matching not found");
  }
  change_res(21, F.dual[21], false, S);

  change_res(22, F.dual[22], true, S);
  if (!assign_match_edges(0, F, S)) {
    throw_error(F.n, h, p, 1, "Perfect matching not found");
  }
  change_res(22, F.dual[22], false, S);
  change_res(24, F.dual[24], false, S);

  // if the face 11 is chosen,
  // it will work with all faces from patch_2
  change_res(11, F.dual[11], true, S);
  for (int i = 0; i < 6; i++) {
    p_2 = pent_patch_2[i];
    change_res(p_2, F.dual[p_2], true, S);
    if (!assign_match_edges(0, F, S)) {
      throw_error(F.n, h, p, 1, "Perfect matching not found");
    }
    change_res(p_2, F.dual[p_2], false, S);
  }
  // all faces but 0 will work from patch_1
  for (int i = 1; i < 5; i++) {
    p_2 = pent_patch_1[i];
    change_res(p_2, F.dual[p_2], true, S);
    if (!assign_match_edges(0, F, S)) {
      throw_error(F.n, h, p, 1, "Perfect matching not found");
    }
    change_res(p_2, F.dual[p_2], false, S);
  }
  change_res(0, F.dual[0], true, S);
  if (assign_match_edges(0, F, S)) {
    throw_error(F.n, h, p, 1, "Perfect matching found");
  }
  change_res(0, F.dual[0], false, S);
  change_res(11, F.dual[11], false, S);

  cout << "Unit tests successful" << endl;
}
