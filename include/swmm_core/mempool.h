//-----------------------------------------------------------------------------
//  mempool.h
//
//  Header for mempool.c
//
//  The type alloc_handle_t provides an opaque reference to the
//  alloc pool - only the alloc routines know its structure.
//-----------------------------------------------------------------------------
#ifndef MEMPOOL_H
#define MEMPOOL_H




struct alloc_handle_t
{
	long  dummy;
};
typedef struct alloc_handle_t alloc_handle_t;

typedef struct Project Project;

alloc_handle_t *AllocInit(Project*);
char           *Alloc(Project*,long);
alloc_handle_t *AllocSetPool(Project*, alloc_handle_t *);
void            AllocReset(Project*);
void            AllocFreePool(Project*);

#endif
