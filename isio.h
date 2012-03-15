/*
   March 2012 - Adapted for use with samtools by Kevin Lewis (Wellcome Trust
   Sanger Institute) from code supplied under a BSD license by the iRODS
   project (www.irods.org).
*/

/* iRODS unix stdio i/o functions; emulate a subset of stdio.h
 * calls */

#ifndef IRODS_IO_H
#define IRODS_IO_H

#include <stdio.h>

typedef struct _isio_file {
	long base_offset;    /* how far base is from start of file */
	char *base;          /* actual data cache */
	int bufferSize;      /* size of cache (actual space malloced) */
	char *ptr;           /* current read position in cache */
	int count;           /* how many unread bytes remain in cache */
	int dirty;          /* set to true to indicate buffer contents differ from iRODS storage */
	long l1descInx;	/* returned by rcDataObjCreate() */
} ISIO_FILE;

FILE *irodsfopen(char *filename, char *modes);
size_t irodsfread(void *buffer, size_t itemsize, size_t nitems, FILE *fi_stream);
size_t irodsfwrite(void *buffer, size_t itemsize, size_t nitems, FILE *fi_stream);
int irodsfclose(FILE *fi_stream);
long irodsftell(FILE *fi_stream);
int irodsfseek(FILE *fi_stream, long offset, int whence);
int irodsfflush(FILE *fi_stream);
int irodsfputc(int inchar, FILE *fi_stream);
int irodsfgetc(FILE *fi_stream);
void irodsexit(int exitValue);

ISIO_FILE *isioFileOpen(const char *filename, char *modes);
int isioFillBuffer(ISIO_FILE *ifp);
int isioFileRead(ISIO_FILE *ifp, void *buffer, size_t maxToRead);
int isioFileWrite(ISIO_FILE *ifp, void *buffer, size_t countToWrite);
int isioFileClose(ISIO_FILE *ifp);
int isioFileTell(ISIO_FILE *ifp);
int isioFileSeek(ISIO_FILE *ifp, long offset, int whence);
int isioFlush(ISIO_FILE *ifp);
int isioFilePutc(int inchar, ISIO_FILE *ifp);
int isioFileGetc(ISIO_FILE *ifp);

#endif /* IRODS_IO_H */
