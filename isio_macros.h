/*
   March 2012 - Adapted for use with samtools by Kevin Lewis (Wellcome Trust
   Sanger Institute) from code supplied under a BSD license by the iRODS
   project (www.irods.org).
*/

/* iRODS unix stdio i/o functions; emulate a subset of stdio.h
 * calls */

#ifndef ISIO_MACROS_H
#define ISIO_MACROS_H

#include <stdio.h>

#include "isio.h"

#define fopen(A,B) irodsfopen(A,B)
#define fread(A,B,C,D) irodsfread(A,B,C,D)
#define fclose(A) irodsfclose(A)
#define exit(A) irodsexit(A)
#define fwrite(A,B,C,D) irodsfwrite(A,B,C,D)
#define fseek(A,B,C) irodsfseek(A,B,C)
#define ftell(A) irodsftell(A)
#define fflush(A) irodsfflush(A)
#define fputc(A, B) irodsfputc(A,B)
#define fgetc(A) irodsfgetc(A)

#endif /* IRODS_IO_H */
