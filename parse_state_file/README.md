# `state_to_nibble`

Convert `IQtree` state files and input sequence files to nibble encoded sequences.

## Synopsis

`state_to_nibble <leaf.fa> <nodes.state> <out.prefix>`

`<leaf.fa>` A fasta file of the aligned sequences used as input for `IQtree`

`<nodes.state>` Ancestral state reconstructions output by `IQtree` when using the `-asr` option.

Both the ancestral and extant sequences will be converted to the nibble format and encoded in
`<prefix>.nb_states`. In addition the the maximum likelihood associated with a given location
and node will be quantised to 256 levels and encoded within the `<prefix>.nb_lkhood` file in
a similar manner.

## Compilation

Modify (if desired) and run `compile.sh` to compile the utility. Note
that the script assumes the presence of `bash` at `/bin/bash`.

## Bugs?

Note that the program has not been extensively tested and that not all output files will work.

In particular, the program does not accept a series of likelihoods that sum to more than one.
Such rows may be present due to rounding errors and this constraint should probably be moderated.


## Output formats

### `.nb_states`

-------------     ---------------------------------------------------------------------------------
Bytes             Description
-------------     ---------------------------------------------------------------------------------
1-8               'iqstates'

9-12               The length of each sequence. The sequences should be aligned, and hence all 
                   must have the same length.
				   
13-17              The number of nodes.

18-                The extant and ancestral states. Eight bases are encoded within unsigned 32 bit
18+4*l*n           integers. The actual length of this block is:  
                   `4*l*n + (l % 8 == 0 ? 0 : 4)`

18+4+l*n + 1       The node and leaf names as 0 delimited character arrays.
to end
----------------------------------------------------------------------------------------------------

The reason for using integer encoding is simply because this makes it possible to read the data
into `R` using `readBin` without having do do much fancy parsing.

### `.nb_lkhood`

As `nb_states`, except that each base requires one byte.

