//---------------------------------------------------------------------------
//
// merge class
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// (C) 2011 by Argonne National Laboratory.
// See COPYRIGHT in top-level directory.
//
//----------------------------------------------------------------------------

#include "merge.hpp"

//----------------------------------------------------------------------------
//
// constructor
//
// comm: MPI communicator
//
Merge::Merge(MPI_Comm comm) {

  this->comm = comm;

}
//----------------------------------------------------------------------------
//
// computes group number and position within that group for my block
// to participate in the merge communication
//
// group number is 0 to the global number of groups in the current round - 1
// position number is 0 to k value of the current round - 1
//
// inputs:
// r: current round
// kv: vector of k values
// gid: global id of the block
//
// outputs:
// g: group number
// p: position number within the group
//
inline void Merge::GetGrpPos(int r, const int *kv, int gid, 
			     int &g, int &p) {

  int step = 1;

  for (int i = 0; i < r; i++)
    step *= kv[i];

  // the second term in the following expression does not simplify to
  // gid / kv[r]
  // because the division gid / (step * kv[r]) is integer and truncates
  // this is exactly what we want
  g = gid % step + gid / (step * kv[r]) * step;

  p = gid / step % kv[r];

}
//----------------------------------------------------------------------------
//
// gets the global ids of the blocks in my group
//
// inputs:
// kv: vector of k values
// cur_r: current round number (0 to r - 1)
// gid: global id of the block
//
// output:
// partners: global ids of the partners (blocks) in my group, including myself
// grp: global communication group number
//
// returns: true: my gid is the root (maximum) of the group, else false
//
inline bool Merge::GetPartners(const int *kv, int cur_r, 
			       int gid, int *partners, int grp) {

  int step = 1; // gids jump by this much in the current round
  int p; // position of this gid in the group
  int r, k;

  GetGrpPos(cur_r, kv, gid, grp, p);

  for (r = 0; r < cur_r; r++)
    step *= kv[r];

  partners[0] = gid - p * step;
  for (k = 1; k < kv[cur_r]; k++)
    partners[k] = partners[k - 1] + step;

  return(gid == partners[kv[cur_r] - 1]);

}
//----------------------------------------------------------------------------
//
// radix-k merge
//
// it_in: pointers to input items
// hdrs: pointers to input headers
// nb_in: number of input items
// it_out: pointers to output items
// nr: number of rounds
// kv: k vector, radix for each round
// io: pointer to io object
// assign: pointer to assignment object
// merge_func: pointer to merging function
// item_func: pointer to function that creates item
// type_func: pointer to function that creates MPI datatype for item 
//
// side effects: allocates output items and array of pointers to them
//
// returns: number of output items
//
int Merge::MergeBlocks(char **it_in, int **hdrs, int nb_in, char** &it_out,
		       int nr, int *kv, IO *io, Assignment *assign,
		       char *(*merge_func)(char **, int), 
		       char * (*item_func)(int *),
		       MPI_Datatype* (*type_func)(char *)) {

  int partners[MAX_K]; // ranks of the k procs in my group, including myself
  int rank, groupsize; // MPI usual
  int gid; // global id of current item block
  int p; // process rank
  double t[MAX_R]; // time for each round
  vector<int> gids; // distinct gids of root blocks
  int unused; // argument for a function call but we don't use it

  // init
  assert(nr > 0 && nr <= MAX_R); // sanity
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);

  // work with a copy of it_in, don't disturb original
  vector<char *> my_it_in(it_in, it_in + nb_in);
  vector<bool> done(nb_in, false);

  // for all rounds
  for (int r = 0; r < nr; r++){

    t[r] = MPI_Wtime();
    int n_recv = 0; // number of received blocks by root
    gids.clear();
    
    // all my current blocks must participate in a round (send or receive)
    for (int b = 0; b < nb_in; b++) {

      if (!done[b]) { // skip blocks that did not survive to this round
	gid = assign->RoundRobin_lid2gid(b);
	bool root = GetPartners(kv, r, gid, partners, unused);
	if (!root) { // nonroots post sends of headers and MSCs
	  p = assign->RoundRobin_gid2proc(partners[kv[r] - 1]);
	  // root gid (3rd arg) is the group identifier
	  MPI_Datatype *dtype = type_func(my_it_in[b]);
	  io->SendItem(my_it_in[b], hdrs[b], p, partners[kv[r] - 1], dtype);
	  MPI_Type_free(dtype);
	  done[b] = true; // nonroot blocks are done after they have been sent
	}
	else { // root posts receives of headers
	  gids.push_back(partners[kv[r] - 1]);
	  for (int k = 0; k < kv[r] - 1; k++) { // receive the others
	    p = assign->RoundRobin_gid2proc(partners[k]);
	    // root gid (2nd arg) is the group identifier
	    io->StartRecvItem(p, partners[kv[r] - 1]);
	    n_recv++;
	  }
	}
      }

    }

    // finish receiving all items
    char **mhs = new char*[nb_in * kv[r]];
    int *ids = new int[gids.size() * kv[r]]; // gids of the mhs
    io->FinishRecvItems(mhs, ids, item_func, type_func); 

    // root does a merge for each group it owns
    for (int j = 0; j < (int)gids.size(); j++) {
      vector<char *>mh_temp; // temp. copy of items
      int lid = assign->RoundRobin_gid2lid(gids[j]);
      mh_temp.push_back(my_it_in[lid]);
      for (int i = 0; i < n_recv; i++) { // collect items for this group
	if (ids[i] == gids[j])
	  mh_temp.push_back(mhs[i]);
      }
      my_it_in[lid] = merge_func(&mh_temp[0], kv[r]);
      for (int i = 0; i < kv[r]; i++) // cleanup
	delete mh_temp[i];
    }

    delete[] ids;
    delete[] mhs;
    t[r] = MPI_Wtime() - t[r];

  } // for all rounds

  // copy out and cleanup
  int nb_out = gids.size();
  it_out = new char*[nb_out];
  for (int i = 0; i < nb_out; i++) {
    int lid = assign->RoundRobin_gid2lid(gids[i]);
    it_out[i] = my_it_in[lid];
  }

#if 0
  // debug: print round times
  for (int r = 0; r < nr; r++)
    fprintf(stderr, "Merge round %d took %.3lf s\n", r, t[r]);
#endif

  return nb_out;

}
//----------------------------------------------------------------------------
