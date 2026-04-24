#ifndef _IQ_TOOLS_COMMON
#define _IQ_TOOLS_COMMON

/// nibble: an uint32_t holding up to 8 nibbles
//  i: the last position in the sequence (0 based!)
// Shifts nibbles added since the last write to the left
// of the nibble integer.
#define END_NIBBLES(nibble, i) ( (nibble) <<= (4 * (8 - (i) % 8)) )

// As END_NIBBLES, but for 8 bit encoded quality values instead:
#define END_QUAL(qual, i) ( (qual) <<= (8 * (4 - (i) % 4)) )


#endif
