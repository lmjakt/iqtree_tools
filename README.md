# IQTree utilities

A small number of tools for handling data from IQtree runs. In particular to
facilitate the analysis of ancestral state transitions predicted by IQTree.

Currently there are two tools:

1. `parse_state_file/state_to_nibble`
2. `extract_transitions/extract_transitions`

The first of these is used to extract the most likely ancestral states
from IQTree `.state` files and extant sequences into a more compact
nibble representation. The second tool uses the files output by 
`state_to_nibble` along with the `.treefile` produced by IQTree
to count the transitions at each node.

Both tools are in development and may be changed in the future.

For documentation see the `README.md` files and the source code.

