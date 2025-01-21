#include "include.h"

void change_res(const int f_id, const face f_info, bool res, Clar_struct(&S)) {
  const int f_size = f_info.size; // grab size of face
  // if face is now resonant, everything will be shifted by +1,
  // if face is now non-resonant, everything will be shifted by -1
  const int shift = 2 * res - 1;

  // face is now resonant
  if (res) {
    // note that S.num_res_p + S.num_res_h acts as an index
    S.res_f[S.num_res_p + S.num_res_h] = f_id; // add to resonant set
  }

  // update number of resonant faces appropriately
  if (f_size == 5)
    S.num_res_p += shift; // is resonant pentagon
  else
    S.num_res_h += shift; // is resonant hexagon
  // for each of f_id neighbors and vertices
  for (int i = 0; i < f_size; i++) {
    S.assigned_f[f_info.adj_f[i]] += shift;
    S.covered_v[f_info.vertices[i]] += shift;
  }
}

bool face_term_cond_met(int *f_id, int f_count, const Fullerene(&F),
                        Clar_struct(&S), const int h, const int p,
                        const ofstream out_files_ptr[NFILE]) {
  // correct number of pentagons have been assigned
  if (p == S.num_res_p) {
    // correct number of hexagons have been assigned
    if (h == S.num_res_h) {
#if DEBUG_CLAR
      print_vec(S.res_f, "Face assignment complete, resonant faces: ");
#endif
      // assign matching edges
      return true; // terminate
      // need to assign more resonant hexagons, if there are more to check
    } else if (f_count < F.dual_n)
      return false;
    return true; // no more hexagons to consider, terminate
    // recall faces are in sorted order with pentagons at front.
    // if there are more pentagons to consider
  } else if (f_count < 12)
    return false; // do not terminate
  return true;    // no more pentagons to consider, terminate
}

void assign_res_face(int *f_id, int f_count, const Fullerene(&F),
                     Clar_struct(&S), const int h, const int p,
                     const ofstream out_files_ptr[NFILE]) {
  // if the termination condition has been met for the faces, stop
  if (face_term_cond_met(f_id, f_count, F, S, h, p, out_files_ptr))
    return;

#if DEBUG_CLAR
  cout << "Considering face id: " << *f_id << ", f_count: " << f_count << endl;
#endif

  // termination condition was not met, we need to assign more resonant faces
  // check if f_id has already been assigned (and is therefore non-res)
  if (S.assigned_f[*f_id] > 0) {
    // move onto next face
    assign_res_face(f_id + 1, f_count + 1, F, S, h, p, out_files_ptr);
    return;
  }

  // if the face can be assigned as resonant, make it resonant
  if ((S.num_res_p < p && F.dual[*f_id].size == 5) ||
      (S.num_res_h < h && F.dual[*f_id].size == 6)) {
#if DEBUG_CLAR
    cout << "Face: " << *f_id << " is now resonant" << endl;
#endif
    change_res(*f_id, F.dual[*f_id], true, S);
    // proceed to next face
    assign_res_face(f_id + 1, f_count + 1, F, S, h, p, out_files_ptr);
#if DEBUG_CLAR
    cout << "Face: " << *f_id << " is now non-resonant" << endl;
#endif
    change_res(*f_id, F.dual[*f_id], false, S);
  }

  // proceed to next face with f_id as non-resonant
  assign_res_face(f_id + 1, f_count + 1, F, S, h, p, out_files_ptr);
  return;
}

// sort faces so that pentagons come before hexagons and then by ids
bool compare_face(const vector<face>(&dual), const int a, const int b) {
  if (dual[a].size == dual[b].size) {
    return (a < b);
  } else {
    return (dual[a].size < dual[b].size);
  }
}

void anionic_clar_struct_handler(const Fullerene(&F), Clar_struct(&S),
                                 const int h, const int p,
                                 const ofstream out_files_ptr[NFILE]) {
  S.Fix_num_vert(F.n); // resize based on size of graph
  S.Reset_vals();      // reset values

  // vector of face ids, these will be sorted
  vector<int> sorted_f(F.dual_n);
  for (int i = 0; i < sorted_f.size(); i++)
    sorted_f[i] = i;
  // custom sort of face ids
  sort(sorted_f.begin(), sorted_f.end(), [&](int a, int b) { // Sort faces
    return compare_face(F.dual, a, b);
  });
#if DEBUG_CLAR
  print_vec(sorted_f, "Sorted faces: ");
#endif

  assign_res_face(&sorted_f[0], 0, F, S, h, p, out_files_ptr);
}
