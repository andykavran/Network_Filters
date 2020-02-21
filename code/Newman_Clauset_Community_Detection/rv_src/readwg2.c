// Functions to read a network stored in a WG2 file into a NETWORK struct
//
// Mark Newman  14 AUG 06
//
// To use this software, include "network.h" and "readwg2.h" in your program
// and then call the following.
//
// Function calls:
//   int read_network(NETWORK *network, FILE *stream)
//     -- Reads a network from the FILE pointed to by "stream" into the
//        structure "network".  For the format of NETWORK structs see file
//        "network.h".  Returns 0 if read was successful.
//   void free_network(NETWORK *network)
//     -- Destroys a NETWORK struct again, freeing up the memory


// Inclusions

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "network.h"

// Constants

#define LINELENGTH 1000000

// Types

typedef struct line {
  char *str;
  struct line *ptr;
} LINE;

// Globals

LINE *first;
LINE *current;


// Function to read one line from a specified stream.  Return value is
// 1 if an EOF was encountered.  Otherwise 0.

int read_line(FILE *stream, char line[LINELENGTH])
{
  if (fgets(line,LINELENGTH,stream)==NULL) return 1;
  line[strlen(line)-1] = '\0';   // Erase the terminating NEWLINE
  return 0;
}


// Function to read in the whole file into a linked-list buffer, so that we
// can do several passes on it, as required to read the GML format
// efficiently

int fill_buffer(FILE *stream)
{
  int length;
  char line[LINELENGTH];
  LINE *previous;

  if (read_line(stream,line)!=0) {
    first = NULL;                // Indicates empty buffer
    return 1;
  }
  length = strlen(line) + 1;
  first = malloc(sizeof(LINE));
  first->str = malloc(length*sizeof(char));
  strcpy(first->str,line);

  previous = first;
  while (read_line(stream,line)==0) {
    length = strlen(line) + 1;
    previous->ptr = malloc(sizeof(LINE));
    previous = previous->ptr;
    previous->str = malloc(length*sizeof(char));
    strcpy(previous->str,line);
  }
  previous->ptr = NULL;          // Indicates last line

  return 0;
}


// Function to free up the buffer again

void free_buffer()
{
  LINE *thisptr;
  LINE *nextptr;

  thisptr = first;
  while (thisptr!=NULL) {
    nextptr = thisptr->ptr;
    free(thisptr->str);
    free(thisptr);
    thisptr = nextptr;
  }
}


// Function to reset to the start of the buffer again

void reset_buffer()
{
  current = first;
}


// Function to get the next line in the buffer.  Returns 0 if there was
// a line or 1 if we've reached the end of the buffer.

int next_line(char line[LINELENGTH])
{
  if (current==NULL) return 1;
  strcpy(line,current->str);
  current = current->ptr;
  return 0;
}


// Function to count the vertices in the buffer.  Returns number of
// vertices.  Works by simply counting the number of (nonempty) lines in
// the buffer.

int count_vertices()
{
  int count=0;
  char line[LINELENGTH];

  reset_buffer();
  while (next_line(line)==0) if (strlen(line)>0) count++;
  return count;
}


// Function to allocate space for the network structure stored in the
// buffer, determine the parameters (label, degree) of each of the
// vertices, and then read in the adjacency list.

void create_network(NETWORK *network)
{
  int i,j;
  int inc;
  char *ptr;
  char line[LINELENGTH];
  char label1[LINELENGTH];
  char label2[LINELENGTH];

  // Count the vertices

  network->nvertices = count_vertices();

  // Make space for the vertices

  network->vertex = calloc(network->nvertices,sizeof(VERTEX));

  // Go through the file reading the details of each vertex one by one

  reset_buffer();
  for (i=0; i<network->nvertices; i++) {

    // Skip to next nonempty line

    do {
      next_line(line);
    } while (strlen(line)==0);

    // Read in the details of this vertex

    sscanf(line,"%i %s %s %i%n",
	   &network->vertex[i].id,label1,label2,&network->vertex[i].degree,
	   &inc);

    // Store the label

    strcat(label1," ");
    strcat(label1,label2);
    network->vertex[i].label = malloc((strlen(label1)+1)*sizeof(char));
    strcpy(network->vertex[i].label,label1);

    // Now read in the edges

    network->vertex[i].edge = malloc(network->vertex[i].degree*sizeof(EDGE));

    ptr = line + 2;
    for (j=0; j<network->vertex[i].degree; j++) {
      ptr += inc;
      sscanf(ptr," (%i, %lf)%n",
	     &network->vertex[i].edge[j].target,
	     &network->vertex[i].edge[j].weight,
	     &inc);
    }

  }
}


// Function to read a complete network

int read_network(NETWORK *network, FILE *stream)
{
  fill_buffer(stream);
  create_network(network);
  free_buffer();

  return 0;
}


// Function to free the memory used by a network again

void free_network(NETWORK *network)
{
  int i;

  for (i=0; i<network->nvertices; i++) {
    free(network->vertex[i].edge);
    free(network->vertex[i].label);
  }
  free(network->vertex);
}
