#include <zlib.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include "kstring.h"
#include "kvec.h"
#include "kseq.h"
#include "../common/common.h"

// END_NIBBLES, END_QUAL are defined in common.h

// Use kseq's ks_getuntil() to read line by line as in the example at:
// https://attractivechaos.github.io/klib/#Kseq%3A%20stream%20buffer%20and%20FASTA%2FQ%20parser
//
// Note that kseq defines kstring_t, so it is not necessary to include kstring.h
// for this. If it is included then it is also necessary to compile and link kstring.c
//
// I use gzfile as this can be used with both compressed and noncompressed files.
//
// This tool will pick the most likely sequence from a state file produced by
// IQtree. The sequence will use the four bytes of an unsigned integer to encode
// 8 nibbles, and 4 bytes of one integer to encode 4 likelihood values quantised
// to 256 levels (0 -> 255) / 255.
//
// The tool will write this to two files consisting of:
//
// File 1: extension ".nb_states"
// 1. Two 32 bit integers of the value 'iqstates' (using ASCII encoding)
//    iqst = 0x69717374 ates=0x61746573
// 2. A single unsigned 32 bit integer giving the length (l) of each sequence;
// 3. A single unsigned 32 bit integer giving the number of nodes (nodes_n)
// 3. nodes_n * (l / 8) + (l % 8 > 0 ? 1 : 0) 32 bit integers specifying the sequence
//    in nibble encoding.
// File 2: extension ".nb_lkhood"
// 1. Two 32 bit integers of the value 'iqlkhood' (using ASCII encoding)
//    iqlk =0x69716C6B    hood=0x686F6F64
// 2. A single unsigned 32 bit integer giving the length (l) of each sequence;
// 3. A single unsigned 32 bit integer giving the number of nodes (nodes_n)
// 4. nodes_n * (l / 4) + (l % 4 > 0 ? 1 : 0) 32 bit integers specifying the quantized
//    likelihood values of the predicted states.
//
// 
// 
// It is an error for any sequence to not have the same length.

// The rationale for using this encoding is that it can easily be read into R using
// readBin(), and then interrogated using bitwise operations on the resulting values.
// R does not allow bitwise operations on character values and only has 32 bit integer values
// (which are treated as signed, but as long as we don't do arithmetic with these then
// that doesn't matter).

// The IQtree state file format assumed has three sections:
// 1. Comments preceded by "#" characters. It will be assumed these are only found at the
//    beginning of lines
// 2. A single line header (tab delimited). Should contain
//    Node Site State p_A p_C p_G p_T ## total of 7
// 3. Data lines tab delimited. Columns as in the header.

// KSTREAM_INIT(gzFile, gzread, 16384)
KSEQ_INIT(gzFile, gzread)
  
// These were taken from:
// https://www.bioinformatics.org/sms/iupac.html
// (- and . indicate gaps
// R = A,G 0x1 + 0x4 = 0x5
// Y = C,T 0x2 + 0x8 = 0xA
// S = G,C = 0x2 + 0x4 = 0x6
// W = A,T = 0x1 + 0x8 = 0x9
// K = G,T = 0x4 + 0x8 = 0xC
// M = A,C = 0x1 + 0x2 = 0x3
// B = C,G,T = 0x2+0x4+0x8 = 0xE
// D = A,G,T = 0x1+0x4+0x8 = 0xD
// H = A,C,T = 0x1+0x2+0x8 = 0xB
// V = A,C,G = 0x1+0x2+0x4 = 0x7
// N = A,C,G,T = 0x1+0x2+0x4+0x8 = 0xF

// This can be done more efficiently by converting more than one nucleotide at
// a time as in htslib.
uint8_t c_to_nibble[256] = {0};

uint8_t iupac[17] = {'A', 'C', 'G', 'T',
		     'R', 'Y', 'S', 'W',
		     'K', 'M', 'B', 'D',
		     'H', 'V', 'N',
		     '-', '.'};
uint8_t nibble[17] = {0x1, 0x2, 0x4, 0x8,
		      0x5, 0xA, 0x6, 0x9,
		      0xC, 0x3, 0xE, 0xD,
		      0xB, 0x7, 0xF,
		      0x0, 0x0};
void set_nibble()
{
  for(int i=0; i < 17; ++i){
    c_to_nibble[ iupac[i] ] = nibble[i];
    // A -> 0x41, a -> 0x61
    // tolower -> OR with 0x20
    if(i < 16)
      c_to_nibble[ iupac[i] | 0x20 ] = nibble[i];
  }
}

typedef kvec_t(kstring_t) node_list;

void set_kstring(kstring_t *dest, const char* source, size_t l){
  dest->l = 0;
  kputsn(source, l, dest);
}

// source must be a null terminated string
void set_kstring_0(kstring_t *dest, const char* source){
  dest->l = 0;
  kputsn(source, strlen(source), dest);
}

// fname: fasta file of aligned sequences
// state_fd: file to write states to
// lhood_fd: file to write likelihoods
// seq_length: the expected sequence length. All sequences should be of this
//             lenth.
// nodes: a kvec_t of node names. The sequences names will be appended.
int read_write_leaves(const char *fname, FILE *state_fd, FILE *lhood_fd,
			 uint32_t seq_length, node_list *nodes){
  gzFile fp = gzopen(fname, "r");
  assert(fp != 0 && "Failed to open leaf sequence file\n");
  kseq_t *seq = kseq_init(fp); 
  int l = 0;
  kstring_t null_node = {0, 0, NULL};
  int error = 0;
  while((l = kseq_read(seq)) >=0){
    if(l != seq_length){
      error = -1;
      break;
    }
    kv_push( kstring_t, *nodes, null_node );
    set_kstring(&nodes->a[nodes->n-1], seq->name.s, seq->name.l);
    uint32_t nibble = 0;
    uint32_t qual = 0;
    for(size_t i=0; i < seq->seq.l; ++i){
      nibble = (nibble << 4) | c_to_nibble[ (uint8_t)seq->seq.s[i] ];
      qual = (qual << 8) | 0xFF;
      if(((i+1) % 8) == 0)
	assert( fwrite(&nibble, sizeof(uint32_t), 1, state_fd) == 1 );
      if(((i+1) % 4) == 0)
	assert( fwrite(&qual, sizeof(uint32_t), 1, lhood_fd) == 1);
    }
    if((seq->seq.l % 8) > 0){
      END_NIBBLES(nibble, seq->seq.l);
      assert( fwrite(&nibble, sizeof(uint32_t), 1, state_fd) == 1 );
    }
    if((seq->seq.l % 4) > 0){
      END_QUAL(qual, seq->seq.l);
      assert( fwrite(&qual, sizeof(uint32_t), 1, lhood_fd) == 1);
    }
  }
  kseq_destroy(seq);
  gzclose(fp);
  return(error);
}

int main(int argc, char *argv[])
{
  
  if(argc != 4){
    fprintf(stderr, "Usage: %s <leaf.seq> <nodes.state> <out.prefix>\n", argv[0]);
    return(1);
  }
  set_nibble();

  // File objects and associated.
  gzFile nodes_fp;
  kstream_t *nodes_ks;
  kstring_t line = {0, 0, NULL};
  
  nodes_fp = gzopen(argv[2], "r");
  assert(nodes_fp != 0 && "File pointer set to 0\n");
  
  nodes_ks = ks_init(nodes_fp);

  kstring_t states_out = {0, 0, NULL};
  kstring_t lkhood_out = states_out;

  // Set up the output file names
  set_kstring_0(&states_out, argv[3]);
  set_kstring_0(&lkhood_out, argv[3]);
  kputs(".nb_states", &states_out);
  kputs(".nb_lkhood", &lkhood_out);
  
  // Keep track of node identifiers:
  kstring_t last_node = {0, 0, NULL};
  kstring_t node = {0, 0, NULL};
  kstring_t null_node = {0, 0, NULL};
  //  kvec_t(kstring_t) nodes;
  node_list nodes;
  kv_init(nodes);

  // These are 32 bit values for compatability with R. 
  // A potential problem in the future.
  uint32_t seq_l = 0;
  uint32_t seq_i = 0;
  size_t line_no = 0;
  
  uint8_t state = 0;
  float max_p = 0;
  uint32_t nibble = 0;
  uint32_t nibble_lk = 0;

  // it would be clearer to writ this as
  // char state_magic[8] = {'i', 'q', 's', 't', 'a', 't', 'e', 's'};
  uint32_t state_magic[4] = {0x74737169, 0x73657461, 0, 0};
  uint32_t lhood_magic[4] = {0x6B6C7169, 0x646F6F68, 0, 0};
  FILE *state_fd = fopen(states_out.s, "w");
  FILE *lhood_fd = fopen(lkhood_out.s, "w");

  // This will write the magic numbers and two 0s. We will write the
  // length of the sequence and number of nodes after reading the complete
  // file.
  assert( fwrite(&state_magic, 4, sizeof(uint32_t), state_fd) == 4);
  assert( fwrite(&lhood_magic, 4, sizeof(uint32_t), lhood_fd) == 4);

  long seq_l_pos = 2 * sizeof(uint32_t);
  // past_header is used to skip the first line after comments; that
  // line contains column names. We expect that the header should contain
  // "Node    Site    State   p_A     p_C     p_G     p_T"
  // And it would make sense to write a check for this.
  int past_header = 0;
  // offsets will hold pointers to a split string
  // n_fields gives the number of offsets
  int *offsets = 0;
  int n_fields;
  
  while( ks_getuntil(nodes_ks, '\n', &line, 0) >= 0 ){
    ++line_no;
    if(line.l == 0 || line.s[0] == '#')
      continue;
    if(!past_header){
      past_header = 1;
      continue;
    }
    
    // We expect 7 fields in the file:
    offsets = ksplit(&line, '\t', &n_fields);
    
    if(n_fields != 7){
      fprintf(stderr, "Failed to obtain 7 fields from:\n %s\n", line.s);
      return(1);
    }
    set_kstring(&node, line.s, offsets[1]-1);
    if(!last_node.s){
      set_kstring(&last_node, node.s, node.l);
      kv_push(kstring_t, nodes, null_node);
      set_kstring(&nodes.a[nodes.n-1], node.s, node.l);
    }
    if(strcmp(node.s, last_node.s)){
      printf("new node. %s -> %s  seq_i: %u  seq_l: %u  line: %lu\n",
	     last_node.s, node.s, seq_i, seq_l, line_no);
      kv_push(kstring_t, nodes, null_node);
      set_kstring(&nodes.a[nodes.n-1], node.s, node.l);
      if(seq_l == 0)
	seq_l = seq_i;
      if(seq_i != seq_l){ // error
	fprintf(stderr, "Length of sequence %s (%u) != %u  at line: %lu\n",
		node.s, seq_i, seq_l, line_no);
	exit(3);
      }
      // Otherwise write the last sequence position and likelihoods
      // if not already written.
      // Note: the nibbles are arranged in the order: 1,2,3,4,5,6,7,8
      //       within the nibble uint32_t. If one more nibble is added
      //       then is added the resulting nibbles in the integer will be:
      //       2,3,4,5,6,7,8,9
      // If that is the last nibble added, then we will want to shift the 9th
      // nibble to first position (i.e. << 4 * [8 - extra_nibble_no]
      // where extra_nibble_no will be: (seq_i % 8)
      // This expression is easy to get wrong, use a macro instead...
      // 
      if(seq_i % 8){
	END_NIBBLES(nibble, seq_i);
	assert( fwrite(&nibble, sizeof(uint32_t), 1, state_fd) == 1);
      }
      if(seq_i % 4){
	END_QUAL(nibble_lk, seq_i);
	assert( fwrite(&nibble_lk, sizeof(uint32_t), 1, lhood_fd) == 1);
      }
      nibble = 0;
      nibble_lk = 0;
      seq_i = 0;
      set_kstring(&last_node, node.s, node.l);
    }
    ++seq_i;
    
    state = line.s[ offsets[2] ];
    max_p = 0;
    for(size_t i=3; i < 7; ++i){
      float p = atof( line.s + offsets[i] );
      max_p = p > max_p ? p : max_p;
    }
    assert(max_p <= 1 && "Likelihood larger than 1 encountered\n");

    // Update the nibble and the nibble_lk values
    nibble = (nibble << 4) | c_to_nibble[state];
    nibble_lk = (nibble_lk << 8) | ((uint8_t)( max_p * 255.0 ));

    // If the nibble or nibble_lk is complete write to file.
    if((seq_i % 8) == 0){
      assert( fwrite(&nibble, sizeof(uint32_t), 1, state_fd) == 1 );
    }
    if((seq_i % 4) == 0)
      assert( fwrite(&nibble_lk, sizeof(uint32_t), 1, lhood_fd) == 1);
  }
  if(seq_i % 8){
    END_NIBBLES(nibble, seq_i);
    assert( fwrite(&nibble, sizeof(uint32_t), 1, state_fd) == 1);
  }
  if(seq_i % 4){
    END_QUAL(nibble_lk, seq_i);
    assert( fwrite(&nibble_lk, sizeof(uint32_t), 1, lhood_fd) == 1);
  }
  
  int leaf_er = read_write_leaves(argv[1], state_fd, lhood_fd,
				  seq_l, &nodes);
  if(leaf_er != 0){
    fprintf(stderr, "Encountered error in reading leaf sequences %d\n", leaf_er);
    exit(leaf_er);
  }
  
  // then write the node names as 0 delimited strings;
  for(size_t i=0; i < nodes.n; ++i){
    assert(fwrite(nodes.a[i].s, 1, nodes.a[i].l+1, state_fd) == nodes.a[i].l+1);
    assert(fwrite(nodes.a[i].s, 1, nodes.a[i].l+1, lhood_fd) == nodes.a[i].l+1);
  }
  // Seek to (almost) the beginning of the file and set the
  // number of states and nodes.
  fseek( state_fd, seq_l_pos, SEEK_SET );
  fseek( lhood_fd, seq_l_pos, SEEK_SET );

  // This was UNDEFINED BEHAVIOUR:
  // fwrite(&nodes.n, sizeof(uint32_t), 1, state_fd);
  // 
  // because nodes.n is size_t (probably 64 bits)
  // but I'm only writing 32 bits
  // This works on Small-endian architectures as the first four
  // bytes are the least significant.
  // Changed to:
  
  uint32_t nodes_n = (uint32_t)nodes.n;
  fwrite(&seq_l, sizeof(uint32_t), 1, state_fd);
  fwrite(&seq_l, sizeof(uint32_t), 1, lhood_fd);
  fwrite(&nodes_n, sizeof(uint32_t), 1, state_fd);
  fwrite(&nodes_n, sizeof(uint32_t), 1, lhood_fd);
  fclose(state_fd);
  fclose(lhood_fd);
}
