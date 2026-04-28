#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include "kstring.h"
#include "knhx.h"
#include "kvec.h"
#include "../common/common.h"

// if defined then print out data structures after parsing
// to confirm correctness.
#define PRINT_STRUCTS

// Read and parse a newick tree.
// Then traverse tree and output transitions as nibble pairs
// (parent -> child) for each non-root node.
// Output these to a simple format consisting of
// 0..15: "nible_transition" (16 bytes)
// 16..23: length of sequences (uint32_t)
// 24..31: number of nodes (uint32_t)
// 32..35: start of nodes..
//
// Read nibble formatted states as output by state_to_nibble
// for details of the file format see ../parse_state_file
//
// Each transition will be encoded in the lower 16 bits of
// a 32 bit integer. This is to make it easy to analyse the
// the data using R. The program may also output estimated
// transition qualities in the future.

// This holds a series of positions in a file and the names
// associated with those positions.
typedef struct state_fh_list {
  FILE *parent;
  FILE *child;
  uint32_t seq_l;
  uint32_t nodes_n;
  size_t *seq_pos;
  char **node_names;
  char *names_buffer;
} fh_list;

typedef struct tree_data {
  int nodes_n;
  knhx1_t *tree;
  uint32_t *tree_to_states;
} tree_data;

// parent child; the tree indices for the parent and the child
// parent_name, child_name; symbolic names. Note that the struct
//                          does not own these and care must be taken
//                          to not use freed memory.
// counts: an array of all possible nibble state transitions.
//          parent state: bits 5-8
//          child state : bits 1-4
// (0x0 = gap -> gap)
// (0xFF = N -> N)
// (0x23 = C -> G), ...
typedef struct transition_counts {
  int parent, child;
  const char *parent_name, *child_name;
  double d;
  uint32_t counts[256];
} transition_counts;
  
typedef kvec_t(transition_counts) tr_counts_kvec;

fh_list open_state_file(const char *fname){
  fh_list fh;
  fh.parent = fopen(fname, "r");
  fh.child = fopen(fname, "r");
  const size_t magic_l = 8;
  const char *magic = "iqstates";
  char *check = malloc(magic_l);
  assert( fread((void*)check, magic_l, 1, fh.parent) == 1);
  if(strncmp(magic, check, magic_l)){
    fprintf(stderr, "missing magic\n");
    exit(2);
  }
  // free check at the end of the function.
  assert( fread(&fh.seq_l, sizeof(uint32_t), 1, fh.parent) == 1);
  assert( fread(&fh.nodes_n, sizeof(uint32_t), 1, fh.parent) == 1);
  fh.seq_pos = malloc(sizeof(size_t) * fh.nodes_n);
  size_t data_start = 16;
  size_t data_unit_length = sizeof(uint32_t) * (fh.seq_l / 8 + (fh.seq_l % 8 > 0 ? 1 : 0));
  size_t names_start = data_start + data_unit_length * fh.nodes_n;
  for(uint32_t i=0; i < fh.nodes_n; ++i)
    fh.seq_pos[i] = data_start + (i * data_unit_length);
  // calculate the length of the nodes section
  fseek(fh.parent, 0, SEEK_END);
  // Consider using fgetpos which returns fpos_t
  // Or to use fstat() to avoid seeking to the end and back.
  //
  long fsize = ftell(fh.parent);
  assert(fsize > names_start && "File is smaller than expected. Aborting");
  size_t names_l = (size_t)fsize - names_start;
  fh.names_buffer = malloc(names_l);
  fh.node_names = malloc(sizeof(char*) * fh.nodes_n);
  fseek(fh.parent, names_start, SEEK_SET);
  assert(fread(fh.names_buffer, 1, names_l, fh.parent) == names_l);
  size_t ni = 0;
  fh.node_names[ni] = fh.names_buffer;
  for(size_t i = 0; i < names_l-1; ++i){
    if(fh.names_buffer[i] == 0 && (ni+1) < fh.nodes_n)
      fh.node_names[++ni] = fh.names_buffer + i + 1;
  }
#ifdef PRINT_STRUCTS
  for(uint32_t i=0; i < fh.nodes_n; ++i)
    fprintf(stderr, "%u : %s\n", i, fh.node_names[i]);
#endif
  free(check);
  return(fh);
}

// map the tree nodes to the state file handles.
// make an index the stupid way.. this is very ugly..
// the index structure should already be allocated..
int make_tree_to_state_index(fh_list *fh, tree_data *tree){
  int unassigned = 0;
  for(int i=0; i < tree->nodes_n; ++i){
    const char *tree_node = tree->tree[i].name;
    tree->tree_to_states[i] = (uint32_t)-1; // max int value also 0xFFFFFFFF
    for(uint32_t j=0; j < fh->nodes_n; ++j){
      if(strcmp(tree_node, fh->node_names[j]) == 0){
	tree->tree_to_states[i] = j;
#ifdef PRINT_STRUCTS
	fprintf(stderr, "Assigning tree node %d to handle position %u (%s)\n",
		i, j, tree_node);
#endif
	break;
      }
    }
    if(tree->tree_to_states[i] == (uint32_t)-1){
      unassigned++;
      fprintf(stderr, "Unable to assign file position to node %d : %s\n", i, tree_node);
    }
  }
  return(unassigned);
}

tree_data parse_tree(const char *fname){
  FILE *nwk_file = fopen(fname, "r");
  if(!nwk_file){
    fprintf(stderr, "Unable to open %s (newick file)\n", fname);
    exit(2);
  }
  kstring_t nwk_str = {0, 0, NULL};
  if(kgetline(&nwk_str, (kgets_func*)&fgets, nwk_file) != 0){
    fprintf(stderr, "Error reading newick file: %s\n", fname);
    exit(3);
  }
  fclose(nwk_file);
  // Parse the tree:
  int nodes_n = 0;
  int nwk_err = 0;
  knhx1_t *tree = kn_parse( nwk_str.s, &nodes_n, &nwk_err );
#ifdef PRINT_STRUCTS
  fprintf(stderr, "obtained a tree with %d nodes\n", nodes_n);
  // print out a table like structure representing the tree:
  for(int i=0; i < nodes_n; ++i)
    fprintf(stderr, "%d\t%d\t%.3f\t%s\n", i, tree[i].parent, tree[i].d, tree[i].name);
#endif
  tree_data td;
  td.nodes_n = nodes_n;
  td.tree_to_states = calloc( nodes_n, sizeof(uint32_t) );
  td.tree = tree;
  return(td);
}

//
tr_counts_kvec count_transitions(fh_list *fh, tree_data *tree, const char *out_file){
  // this gives the number of uint32_t values to read in from the file.
  size_t read_n = (fh->seq_l % 8 > 0 ? 1 : 0) + fh->seq_l / 8;
  uint32_t *child_buffer = malloc(sizeof(uint32_t) * read_n);
  uint32_t *parent_buffer = malloc(sizeof(uint32_t) * read_n);
  tr_counts_kvec tr_counts;
  kv_init(tr_counts);
  transition_counts counts;
  for(int c=0; c < tree->nodes_n; ++c){
    int p = tree->tree[c].parent;
    if(p < 0)
      continue;
    uint32_t fh_p = tree->tree_to_states[p];
    uint32_t fh_c = tree->tree_to_states[c];
    // Check for unassigned nodes:
    assert( fh_p != (uint32_t)-1 && fh_c != (uint32_t)-1 && "Either parent or child lacks a file handle");
    // fseek returns -1 on error, 0 otherwise. We can assert this
    assert( fseek(fh->parent, fh->seq_pos[fh_p], SEEK_SET) == 0 && "Failed to seek to parent location");
    assert( fseek(fh->child,  fh->seq_pos[fh_c], SEEK_SET) == 0 && "Failed to seek to child location");

    assert(fread(parent_buffer, sizeof(uint32_t), read_n, fh->parent) == read_n &&
	   "Failed to read states from parent" );
    assert(fread(child_buffer, sizeof(uint32_t), read_n, fh->child) == read_n &&
	   "Failed to read states from child");
    memset(&counts.counts, 0, sizeof(uint32_t) * 256);
    counts.parent = p;
    counts.child = c;
    counts.parent_name = tree->tree[p].name;
    counts.child_name = tree->tree[c].name;
    counts.d = tree->tree[c].d;
    for(size_t i=0; i < read_n; ++i){
      for(size_t j=0; j < 8 && (i * 8 + j) < fh->seq_l; ++j){
	counts.counts[(((parent_buffer[i] >> (28 - j*4)) << 4) & 0xF0) |
		      ((child_buffer[i] >> (28 - j*4)) & 0xF)]++;
      }
    }
    kv_push(transition_counts, tr_counts, counts);
  }
  free(child_buffer);
  free(parent_buffer);
  return(tr_counts);
}

void print_summary(tr_counts_kvec counts, const char *fname){
  FILE *fh = fopen(fname, "w");
  assert(fh != 0 && "Unable to open summary file for writing\n");
  for(size_t i=0; i < counts.n; ++i){
    transition_counts c = counts.a[i];
    fprintf(fh, "%d\t%d\t%lf\t%s\t%s",
	    c.parent, c.child, c.d, c.parent_name, c.child_name);
    for(size_t j=0; j < 256; ++j)
      fprintf(fh, "\t%u", c.counts[j]);
    fprintf(fh, "\n");
  }
  fclose(fh);
}

int main(int argc, char *argv[])
{
  if(argc != 4){
    fprintf(stderr, "usage: %s <newick_file> <states_file> <out_prefix>\n", argv[0]);
    exit(1);
  }

  // Obtain the tree
  tree_data tree = parse_tree(argv[1]);
  // Obtain information about the state files:
  fh_list fh_state = open_state_file(argv[2]);

  int unassigned = make_tree_to_state_index(&fh_state, &tree);
#ifdef PRINT_STRUCTS
  fprintf(stderr, "Total of %d nodes in tree were not found in the state file\n",
	  unassigned);
#endif
  
  tr_counts_kvec tr_counts = count_transitions(&fh_state, &tree, "test");

  kstring_t sum_file = {0, 0, NULL};
  kputs(argv[3], &sum_file); kputs("_tr_counts.tsv", &sum_file);
  print_summary(tr_counts, sum_file.s);
}
