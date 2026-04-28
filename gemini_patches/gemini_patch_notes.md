
# extract_transitions.c

## Length of string buffer line, (72)

```diff
-  char *check = malloc(8);
+  char *check = malloc(9); // BUGFIX: Increased allocation to 9 bytes to accommodate null terminator
   assert( fread((void*)check, 8, 1, fh.parent) == 1);
+  check[8] = '\0'; // BUGFIX: Added null terminator to prevent buffer over-read during strcmp
```

Is a legitimate bug because of the use of `strcmp` rather than
`strncmp`. `strcmp` should stop the comparison at the end of `"magic"` but it
doesn't seem to be specified and may be implementation dependent.

Addressed by:

1. Defining variable `magic_l` to 8 (the length of the magic).
2. Use of `magic_l` for calls to `malloc`, `fread` and `strncmp`.

Additionally, not freeing `check` afterwards constitutes a memory leak.
That would only be important if the function is called repeatedly.


## Missing bounds check on file size, (line 88)

The program assumes that the file is longer than a size estimated from the
number of nodes and the length of each sequence. But it doesn't check this.

This is a legitimate bug. Gemini also suggests using `long` rather
than `size_t`. Given that `ftell` returns a `long` this is also correct, and
fine if `long` is 64 bits. It seems that the better approach is to use
`fgetpos` which returns `fpos_t`.

Addressed by use of `assert` and an additional variable `long fsize`.

```diff
-  size_t names_l = ftell(fh.parent) - names_start;
+  long fsize = ftell(fh.parent);
+  if(fsize < names_start) { // BUGFIX: Added bounds check to prevent arithmetic underflow if file is too small
+    fprintf(stderr, "File too small, missing names block\n");
+    exit(2);
+  }
+  size_t names_l = fsize - names_start;

```


## Check value of `tree->tree_to_states[]` for `0xFFFFFFFF`

This is a legitimate bug. I'm also impressed by the fact that it can infer
that I used `(uint32_t)-1` as the default value of
`tree->tree_to_states[p]`.

Addressed by an assert statement. I have also included assert statements to the calls to `fseek`
checking for errors.

```diff
@@ -175,3 +181,5 @@ tr_counts_kvec count_transitions(fh_list *fh, tree_data *tree, const char *out_f
       continue;
     uint32_t fh_p = tree->tree_to_states[p];
     uint32_t fh_c = tree->tree_to_states[c];
+    if(fh_p == (uint32_t)-1 || fh_c == (uint32_t)-1) // BUGFIX: Safely skip nodes unassigned in the states file to prevent segfaults
+      continue;
```

## Changed shift multiplier in extraction

Gemini suggests, `shift multiplier from j*8 to j*4 to prevent negative
shift`. This is actually a serious bug; the multiplier should be 4 as we use 4
bits for nibble encoding. I used 8 bits for quality values. Here the maximum
value of `j*8` will be `7*8` which is far too big.

This bug actually matters as it may result in incorrect output. I'd expect that the
result would be to count half the values twice, with a small number of incorrect
values. But it could also result in a segmentation fault.

```diff
     counts.d = tree->tree[c].d;
     for(size_t i=0; i < read_n; ++i){
       for(size_t j=0; j < 8 && (i * 8 + j) < fh->seq_l; ++j){
-	counts.counts[(((parent_buffer[i] >> (28 - j*8)) << 4) & 0xF0) |
-		      ((child_buffer[i] >> (28 - j*8)) & 0xF)]++;
+	// BUGFIX: Changed shift multiplier from j*8 to j*4 to prevent negative shift amounts (undefined behavior)
+	counts.counts[(((parent_buffer[i] >> (28 - j*4)) << 4) & 0xF0) |
+		      ((child_buffer[i] >> (28 - j*4)) & 0xF)]++;
       }
     }
```

# state_to_nibble.c

## Additional shift at end of sequence

This looks like quite a subtle bug that's quite difficult to find. I'm rather
impressed with the fact that this was identified. The reason it matters is
because the reader expects that the lower positions are found in the more
significant bits. This bug will result in some bases from the prior octet
being repeated and some being missed from the end when the length of the
sequence is not an even multiple of 8.

Fixed by use of macros `END_NIBBLES` and `END_QUAL`.

```diff
-    if((seq->seq.l % 8) > 0)
-      assert( fwrite(&nibble, sizeof(uint32_t), 1, state_fd) == 1 );
-    if((seq->seq.l % 4) > 0)
-      assert( fwrite(&qual, sizeof(uint32_t), 1, lhood_fd) == 1);
+    if((seq->seq.l % 8) > 0) {
+      uint32_t p_nibble = nibble << ((8 - (seq->seq.l % 8)) * 4); // BUGFIX: Shift incomplete sequence left so it aligns with the most significant bits when parsed
+      assert( fwrite(&p_nibble, sizeof(uint32_t), 1, state_fd) == 1 );
+    }
+    if((seq->seq.l % 4) > 0) {
+      uint32_t p_qual = qual << ((4 - (seq->seq.l % 4)) * 8); // BUGFIX: Shift incomplete likelihoods left to match alignment
+      assert( fwrite(&p_qual, sizeof(uint32_t), 1, lhood_fd) == 1);
+    }
```

## Nibble int updated before checking if node is new

This also looks like a real bug. But working out what the fix does exactly is difficult because it
has been split into two different parts; one which adds code and one which removes it.

Instead of accepting the suggestion I have changed the flow of the section to simplify it. It now
does:

For each data line:

1. Split line by tabs.
2. Exit if number of fields is not 7.
3. Set the value of the `node`.
4. If this is the first node, set the value of `last_node` to `node`.
5. If `last_node != node`:
   a. Append `node` to `nodes`
   b. If `seq_l` has not been set, set it to `seq_i`
   c. Exit if `seq_l != seq_i` (all sequences should be the same length)
   d. If `nibble` or `nibble_lk` are partially set, then shift the last
      bits to the left and write to file.
   e. Set `nibble` and `nibble_lk` to 0. This should not be necessary.
   f. Set `seq_i` to 0.
   g. Set `last_node` to `node`.
6. Increment the length counter `seq_i`.
7. Obtain the state and likelihood associated with it.
8. Update `nibble` and `nibble_lk`.
9. If `nibble` or `nibble_lk` is filled then write to file.

At the end, again, check if there is a partial nibble to read or write.

This avoid incrementing `seq_i` prior to the check for a new node and thus avoids
adjusting it in the conditional block.

(The following diff is Gemini's suggestion).

```diff
       set_kstring(&last_node, node.s, node.l);
       kv_push(kstring_t, nodes, null_node);
       set_kstring(&nodes.a[nodes.n-1], node.s, node.l);
+    } else if(strcmp(node.s, last_node.s)){ // BUGFIX: Node transition is checked BEFORE the current line's sequence data is added to the accumulator
+      printf("new node. %s -> %s  seq_i: %u  seq_l: %u  line: %lu\n",
+	     last_node.s, node.s, seq_i-1, seq_l, line_no);
+      nodes_n++;
+      if(seq_l == 0)
+	seq_l = seq_i - 1;
+      if((seq_i - 1) != seq_l){ // error
+	fprintf(stderr, "Length of sequence %s (%u) != %u  at line: %lu\n",
+		last_node.s, seq_i-1, seq_l, line_no);
+	exit(3);
+      }
+      if((seq_i - 1) % 8) {
+	uint32_t p_nibble = nibble << ((8 - ((seq_i - 1) % 8)) * 4); // BUGFIX: Shift partial state sequences left for correct bitwise extraction alignment
+	assert( fwrite(&p_nibble, sizeof(uint32_t), 1, state_fd) == 1);
+      }
+      if((seq_i - 1) % 4) {
+	uint32_t p_lk = nibble_lk << ((4 - ((seq_i - 1) % 4)) * 8);
+	assert( fwrite(&p_lk, sizeof(uint32_t), 1, lhood_fd) == 1);
+      }
+      nibble = 0;
+      nibble_lk = 0;
+      seq_i = 1;
+      set_kstring(&last_node, node.s, node.l);
+      kv_push(kstring_t, nodes, null_node);
+      set_kstring(&nodes.a[nodes.n-1], node.s, node.l);
     }
```

## Division of likelihood used instead of multiplication

Correctly identified.

Corrected as part of addressing the previous issue.

```diff
     nibble = (nibble << 4) | c_to_nibble[state];
-    nibble_lk = (nibble_lk << 8) | ((uint8_t)( max_p / 255 ));
+    nibble_lk = (nibble_lk << 8) | ((uint8_t)( max_p * 255.0 )); // BUGFIX: Changed division to multiplication to prevent likelihood truncation to 0
```

This diff part continues with the removal of code that was replaced by the previous
issue. It has been addressed as described above.

```diff
     if((seq_i % 4) == 0)
       assert( fwrite(&nibble_lk, sizeof(uint32_t), 1, lhood_fd) == 1);
-    
-    if(strcmp(node.s, last_node.s)){
-      printf("new node. %s -> %s  seq_i: %u  seq_l: %u  line: %lu\n",
-	     last_node.s, node.s, seq_i, seq_l, line_no);
-      nodes_n++;
-      kv_push(kstring_t, nodes, null_node);
-      set_kstring(&nodes.a[nodes.n-1], node.s, node.l);
-      // decrement seq_i, because we have gone past the last entry of
-      // the last node.
-      seq_i--;
-      if(seq_l == 0)
-	seq_l = seq_i;
-      if(seq_i != seq_l){ // error
-	fprintf(stderr, "Length of sequence %s (%u) != %u  at line: %lu\n",
-		node.s, seq_i, seq_l, line_no);
-	exit(3);
-      }
-      // Otherwise write the last sequence position and likelihoods
-      // if not already written.
-      // DO CONSIDER CHECKING THE RETURN VALUES!
-      if(seq_i % 8)
-	assert( fwrite(&nibble, sizeof(uint32_t), 1, state_fd) == 1);
-      if(seq_i % 4)
-	assert( fwrite(&nibble_lk, sizeof(uint32_t), 1, lhood_fd) == 1);
-      nibble = 0;
-      nibble_lk = 0;
-      seq_i = 1;
-      set_kstring(&last_node, node.s, node.l);
-      continue;
-    }
   }
-  if(seq_i % 8)
-    assert( fwrite(&nibble, sizeof(uint32_t), 1, state_fd) == 1);
-  if(seq_i % 4)
-    assert( fwrite(&nibble_lk, sizeof(uint32_t), 1, lhood_fd) == 1);
+  if(seq_i % 8) {
+    uint32_t p_nibble = nibble << ((8 - (seq_i % 8)) * 4); // BUGFIX: Pad remaining states for the final sequence
+    assert( fwrite(&p_nibble, sizeof(uint32_t), 1, state_fd) == 1);
+  }
+  if(seq_i % 4) {
+    uint32_t p_lk = nibble_lk << ((4 - (seq_i % 4)) * 8);
+    assert( fwrite(&p_lk, sizeof(uint32_t), 1, lhood_fd) == 1);
+  }
 
   int leaf_er = read_write_leaves(argv[1], state_fd, lhood_fd,
 				  seq_l, &nodes);

```

## Implicit conversion of `size_t` to `uint32_t`

Not identified by Gemini, but noticed due to reviewing the code as
part of the assessing the bug reports.

I had:

```c
fwrite(&nodes.n, sizeof(uint32_t), 1, state_fd);
```

`nodes.n` is `size_t`, but I'm only writing the first 32 bits of it. That works
for small-endian systems, but would fail on a big-endian system.

