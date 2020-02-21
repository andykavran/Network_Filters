// Header file for the parser for WG2 files
//
// Mark Newman  11 AUG 06

#ifndef _READWG2_H
#define _READWG2_H

#include <stdio.h>
#include "network.h"

int read_network(NETWORK *network, FILE *stream);
void free_network(NETWORK *network);

#endif
