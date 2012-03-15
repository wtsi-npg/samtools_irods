/* 
   irods standard i/o emulation library  (initial version)
   
   This package converts Unix standard I/O calls (stdio.h: fopen,
   fread, fwrite, etc) to the equivalent irods calls.  To convert an
   application you just need to add an include statement (isio.h) and
   recompile and relink with irods libraries.  fopen calls with input
   filenames that include the IRODS_PREFIX string ('irods:') will be
   opened and handled as iRODS files; without the prefix this package
   will call the fopen family to handle them as regular local files.

   Like the fopen family, this library does some caching to avoid
   small I/O (network) calls, greatly improving performance.

   The user callable functions are defined in the isio.h and are of
   the form irodsNAME, such as irodsfopen.  Internal function names
   begin with 'isio'.

   The irods environment is assummed.  That is, like i-commands, this
   library needs to read the user's .irodsEnv and authentication files
   to be able to connect to an iRODS server.

   See the ../include/isio.h, ../test*.c files, and ../Makefile for
   more information.

   March 2012 - Adapted for use with samtools by Kevin Lewis (Wellcome Trust
   Sanger Institute) from code supplied under a BSD license by the iRODS
   project (www.irods.org).

 */

#include <stdio.h>
#include <signal.h>

#define NDEBUG 1
#include <assert.h>

#include "rodsClient.h"
#include "dataObjRead.h"

#include "isio.h"

#define IRODS_PREFIX "irods:"
#define ISIO_MAX_OPEN_FILES 20

/* The following two numberic values are also used by the
   test script (modified to be smaller as a test) */
#define ISIO_INITIAL_BUF_SIZE  65536
# define ISIO_MAX_BUF_SIZE    2097152

#define ERRMSSZ 128

int debug=0;

/*
Static function prototypes
*/
static int isioSetup(void);
static ISIO_FILE *new_isio_file(long l1descInx);
static int destroy_isio_file(ISIO_FILE *ifp);
static FILE *add_to_map(ISIO_FILE *ifp);
static ISIO_FILE *is_in_map(int idx);
static int offset_in_current_buffer(ISIO_FILE *ifp, long offset, int whence, long *adj);
static char *copy_filename(const char *filename);

/* set up communcation with iRODS */
static int setupFlag=0;
char localZone[100]="";
rcComm_t *Comm;
rodsEnv myRodsEnv;

/*
isioSetup:
Subsequent calls after successful setup just return 0 (no error)
*/
static int isioSetup(void)
{
	int status;
	rErrMsg_t errMsg;
	char *mySubName;
	char *myName;

	if (debug) printf("isioSetup\n");

	if(setupFlag != 0)
		return(0);

	if((status = getRodsEnv (&myRodsEnv)) < 0) {
		rodsLogError(LOG_ERROR, status, "isioSetup: getRodsEnv error.");
	}
	Comm = rcConnect (myRodsEnv.rodsHost, myRodsEnv.rodsPort, myRodsEnv.rodsUserName, myRodsEnv.rodsZone, 0, &errMsg);

	if(Comm == NULL) {
		myName = rodsErrorName(errMsg.status, &mySubName);
		rodsLog(LOG_ERROR, "rcConnect failure %s (%s) (%d) %s",
				myName,
				mySubName,
				errMsg.status,
				errMsg.msg);
		status = errMsg.status;
		return(status);
	}

	/*
	revert SIGPIPE handler to default.
	Note: this is done to cause sensible behaviour with downstream pipes (head, more,...), but this also means that
	any SIGPIPE signals from the server side will also be fatal. (See _rcConnect() in lib/core/src/rcConnect.c in
	the iRODS source code)
	*/
	signal(SIGPIPE, SIG_DFL);

	if((status = clientLogin(Comm)) == 0) {
		setupFlag=1;
	}

	return(status);
}


/* _cacheInfo handling routines */

static ISIO_FILE *_cacheInfo[ISIO_MAX_OPEN_FILES] = { NULL };

static ISIO_FILE *new_isio_file(long l1descInx)
{
	ISIO_FILE *ret = NULL;

	if((ret=calloc(1, sizeof(ISIO_FILE))) != NULL) {
		ret->base_offset = 0L;
		if((ret->base=malloc(sizeof(char) * ISIO_INITIAL_BUF_SIZE)) == NULL) {
			fprintf(stderr, "new_isio_file: failed to allocate base buffer of size %d\n", ISIO_INITIAL_BUF_SIZE);
			free(ret);
			return(NULL);
		}
		ret->bufferSize = sizeof(char) * ISIO_INITIAL_BUF_SIZE;
		ret->ptr = ret->base;
		ret->count = 0;
		ret->dirty = 0;
		ret->l1descInx = l1descInx;
	}

	assert((ret->base != NULL) && (ret->ptr >= ret->base) && ((ret->ptr - ret->base) + ret->count <= ret->bufferSize));

	return(ret);
}

static int destroy_isio_file(ISIO_FILE *ifp)
{
	if(ifp == NULL)
		return(-1);

	if(ifp->base != NULL)
		free(ifp->base);

	free(ifp);

	return(0);
}

static FILE *add_to_map(ISIO_FILE *ifp)
{
	int i;

	if(ifp == NULL) {
		return(NULL);
	}

	for(i=0; i<ISIO_MAX_OPEN_FILES; i++) {
		if(_cacheInfo[i] == NULL) {
			_cacheInfo[i] = ifp;
			return((FILE *)((unsigned long)i+1));	/* don't use index 0 (would look like NULL) */
		}
	}

	fprintf(stderr, "Too many open files in add_to_map\n");

	return(NULL);
}

static ISIO_FILE *is_in_map(int idx)
{
	idx--;	/* initially 1-based, so adjust to array index */
	if((idx >= 0) && (idx < ISIO_MAX_OPEN_FILES))
		return(_cacheInfo[idx]);
	else
		return(NULL);
}

/*
isioFileOpen:
Return: if error NULL, else pointer to ISIO_FILE

Note: although parseRodsPathStr() treats "filename" as const, the declaration in 
lib/core/include/rodsPath.h doesn't use the "const" keyword.  Since the iRODS code
is unwilling to guarantee no alteration to the contents of filename (and the
compiler complains), I've added the inefficient strdup()/free().
*/
ISIO_FILE *isioFileOpen(const char *filename, char *modes)
{
	int status;
	dataObjInp_t dataObjInp;
	char *fn_copy = NULL;

	if(filename == NULL) return(NULL);

	if(debug) printf("isioFileOpen: %s\n", filename);

	if(isioSetup()) return(NULL);

	memset (&dataObjInp, 0, sizeof (dataObjInp));

	if((fn_copy=copy_filename(filename)) == NULL) return(NULL);
	status = parseRodsPathStr (fn_copy , &myRodsEnv, dataObjInp.objPath);
	free(fn_copy);

	if(status < 0) {
		rodsLogError (LOG_ERROR, status, "isioFileOpen");
		return(NULL);
	}

	/* Set the openFlags based on the input mode (incomplete
		* currently) // */
	dataObjInp.openFlags = O_RDONLY;
	if(strncmp(modes,"w",1)==0) {
		dataObjInp.openFlags = O_WRONLY;
		/* Need to handle other cases sometime */
	}
	if(strncmp(modes,"r+",2)==0) {
		dataObjInp.openFlags = O_RDWR;
	}

	status = rcDataObjOpen (Comm, &dataObjInp);

	if (status==CAT_NO_ROWS_FOUND && dataObjInp.openFlags == O_WRONLY) {
		status = rcDataObjCreate(Comm, &dataObjInp);
	}

	if(status < 0) {
		char errms[ERRMSSZ];

		snprintf(errms, ERRMSSZ, "isioFileOpen(%s, %s)", filename, modes);
		rodsLogError (LOG_NOTICE, status, errms);
		return(NULL);
	}

	return(new_isio_file(status));
}

/*
copy_filename:
	strdup() with some error checking
	Note: the caller is responsible for freeing the allocated memory
*/
static char *copy_filename(const char *filename)
{
	char *fn_copy = NULL;

	if((fn_copy = strdup(filename)) == NULL) {
		rodsLogError(LOG_ERROR, SYS_MALLOC_ERR, "isioFileOpen/copy_filename");
		return(NULL);
	}

	return(fn_copy);
}

FILE *irodsfopen(char *filename, char *modes)
{
	int len;

	if (debug) printf("irodsfopen: %s\n", filename);

	len = strlen(IRODS_PREFIX);
	if(strncmp(filename,IRODS_PREFIX,len)==0) {
		return(add_to_map(isioFileOpen(filename+len, modes)));
	}
	else {
		return(fopen(filename, modes));
	}
}

/*
isioFillBuffer:
Fill ifp cache buffer from current position in file. Adjust ifp->base_offset to new
position in underlying file.
Returns: 0 for success, <0 for error
*/
int isioFillBuffer(ISIO_FILE *ifp)
{
	int status;
	openedDataObjInp_t dataObjReadInp;
	bytesBuf_t dataObjReadOutBBuf;

	if (debug) printf("isioFillBuffer\n");

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	if(ifp->base == NULL)
		return(-1);

	if(ifp->bufferSize <= 0)
		return(-2);

	ifp->base_offset = isioFileTell(ifp);	/* move beginning of cache buffer window to new position */

	dataObjReadOutBBuf.buf = ifp->base;
	dataObjReadOutBBuf.len = ifp->bufferSize;

	memset(&dataObjReadInp, 0, sizeof (dataObjReadInp));

	dataObjReadInp.l1descInx = ifp->l1descInx;
	dataObjReadInp.len = ifp->bufferSize;

	status = rcDataObjRead(Comm, &dataObjReadInp, &dataObjReadOutBBuf);

	if(debug) printf("isioFillBuffer rcDataObjRead stat: %d\n", status);
	if(status < 0) return(status);

/*	ifp->base_offset += status; */
	ifp->ptr = ifp->base;
	ifp->count = status;

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	return(0);
}

/*
isioFileRead:
	Transfer maxToRead bytes from current position in cacheInfo buffer to user buffer.

	If the request cannot be satisfied from the contents of the current cache buffer,
	refill the cache buffer after copying any needed cached data, reallocating the
	buffer if necessary.

	If the read request size (minus any bytes already in the cache) is greater than the
	maximum cacheInfo buffer size of ISIO_MAX_BUF_SIZE, the data should be read directly
	into the user-supplied buffer, and the current cache buffer marked as empty.

	Parameters:
		ifp - data structure for managing previously opened connection to iRODS
		buffer - copy the data to this user-supplied buffer
		maxToRead - the number of bytes to be read; should be <= size of buffer

	Returns:
		number of bytes read, or < 0 to indicate error
*/
int isioFileRead(ISIO_FILE *ifp, void *buffer, size_t maxToRead)
{
	int status;
	int reqSize;	/* number of requested bytes still unread (there may be multiple reads) */
	char *destPtr;	/* destination for data transfers for a given read (within user-supplied buffer) */
	int read_count = 0;	/* number of bytes read, returned to caller */
	int second_read_count = 0;
	int newBufSize;
	int usingUsersBuffer = 0;
	char *saved_base = NULL;
	int saved_bufferSize = 0;

	if (debug) printf("isioFileRead\n");

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	if(maxToRead == 0)
		return(0);

	/* If the buffer had been used for writing, flush it */
	if((status = isioFlush(ifp)) < 0)
		return(status);

	reqSize = maxToRead;
	destPtr = buffer;
	if(ifp->count > 0) {	/* there is something unread in the cache */
		if(ifp->count >= reqSize) {	/* there is enough in the cache for this request */
			memcpy(destPtr, ifp->ptr, reqSize);
			ifp->ptr += reqSize;
			ifp->count -= reqSize;

			return(maxToRead);	/* read request satisfied by cache contents, nothing more to do here */
		}
		else {	/* not enough data in cache for the request - copy what's available and fall through to buffer refill */
			memcpy(destPtr, ifp->ptr, ifp->count);

			reqSize -=  ifp->count;
			read_count = ifp->count;
			destPtr += ifp->count;
			/* mark cache buffer as empty */
			ifp->ptr = ifp->base;
			ifp->count = 0;
		}
	}

	/******************************************************************************
	* At this point, not all data will have been transferred to user buffer because
	* requirements could not be met with cached data. Refill cache buffer and copy
	* remaining requested data to user buffer.
	******************************************************************************/
	newBufSize=(2*reqSize)+8;

	if(ifp->base == NULL)
		ifp->bufferSize = 0;

	/* if necessary, set up a new ifp buffer for the read */
	if(newBufSize > ifp->bufferSize) {
		if(newBufSize<=ISIO_MAX_BUF_SIZE) {
			/* allocate larger cache buffer */
			if((ifp->base=realloc(ifp->base, newBufSize)) == NULL) {
				fprintf(stderr,"Memory Allocation error\n");
				return(-99);
			}

			ifp->bufferSize = newBufSize;

			/* mark cache buffer as empty */
			ifp->ptr = ifp->base;
			ifp->count = 0;
		}
		else {
			/* Read directly into user's buffer. */
			usingUsersBuffer = 1;

			/* Save cache buffer attributes for restoration after read */
			saved_base = ifp->base;
			saved_bufferSize = ifp->bufferSize;

			/* point to user buffer for read */
			ifp->base=destPtr;
			ifp->ptr = ifp->base;
			ifp->bufferSize = reqSize;
		}
	}
	else {
		/* current cache buffer is big enough. Use it, but mark it as empty */
		ifp->ptr = ifp->base;
		ifp->count = 0;

	}

	if((status = isioFillBuffer(ifp)) < 0)
		return(status);

	second_read_count = (ifp->count < reqSize? ifp->count: reqSize);
	read_count += second_read_count;

	if(usingUsersBuffer) {
		/* restore saved cache buffer */
		ifp->base = saved_base;
		ifp->bufferSize = saved_bufferSize;

		/* mark cache buffer as empty */
		ifp->ptr = ifp->base;
		ifp->count = 0;

		if(debug) printf("isioFileRead return1: %d\n", read_count);
	}
	else {
		/* do the copy to the user's buffer */
		if(second_read_count > 0) {
			memcpy(destPtr, ifp->ptr, second_read_count);
			ifp->ptr += second_read_count;
			ifp->count -= second_read_count;
		}

		if(debug) printf("isioFileRead return2: %d\n", read_count);
	}

	return(read_count);
}

size_t irodsfread(void *buffer, size_t itemsize, size_t nitems, FILE *fi_stream)
{
	int i;
	ISIO_FILE *ifp;

	i = (int)(unsigned long)fi_stream;
	if (debug) printf("isiofread: %d\n", i);

	if((ifp=is_in_map(i)) != NULL) {
		return(isioFileRead(ifp, buffer, itemsize*nitems)/itemsize);
	}
	else {
		return(fread(buffer, itemsize, nitems, fi_stream));
	}
}

/*
isioFileWrite:
	Transfer countToWrite bytes from user buffer to current position in cacheInfo buffer,
	flushing to iRODS as needed.

	If the size of the data to be written is more than the remaining space in the current
	cache buffer, refill the cache buffer after flushing the cache buffer to iRODS, reallocating
	the buffer if necessary.

	If the write request size (minus any bytes which fit in the current cache buffer) is greater
	than the maximum cacheInfo buffer size of ISIO_MAX_BUF_SIZE, the data should be read directly
	from the user-supplied buffer, and the current cache buffer marked as empty.

	Parameters:
		ifp - data structure for managing previously opened connection to iRODS
		buffer - copy the data from this user-supplied buffer
		maxToRead - the number of bytes to be written; should be <= size of buffer

	Returns:
		number of bytes read, or < 0 to indicate error
*/
int isioFileWrite(ISIO_FILE *ifp, void *buffer, size_t countToWrite)
{
	int wrtSize;	/* number of bytes to be written in a given write (there may be multiple writes) */
	char *srcPtr;	/* source for data transfers for a given write (within user-supplied buffer) */
	int write_count = 0;	/* number of bytes read, returned to caller */
	int status;
	int spaceInBuffer;
	int newBufSize;
	openedDataObjInp_t dataObjWriteInp;
	bytesBuf_t dataObjWriteOutBBuf;

	if(debug) printf("isioFileWrite\n");

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	if(countToWrite == 0)
		return(0);

	ifp->dirty = 1;

	srcPtr = buffer;
	spaceInBuffer = ifp->bufferSize - (ifp->ptr - ifp->base);
                   
	if (debug) printf("isioFileWrite: spaceInBuffer %d\n", spaceInBuffer);

	wrtSize = (countToWrite <= spaceInBuffer)? countToWrite: spaceInBuffer; /* min(countToWrite, spaceInBuffer) */

	/* write out anything that fits in the cache */
	if (debug) printf("isioFileWrite: caching 1 %p %u\n", ifp->ptr, (unsigned)countToWrite);
	memcpy(ifp->ptr, srcPtr, wrtSize);
	srcPtr += wrtSize;
	ifp->ptr += wrtSize;
	ifp->count -= wrtSize;
	if(ifp->count < 0)
		ifp->count = 0;
	
	write_count = wrtSize;

	if(wrtSize == countToWrite) {
		return(write_count);	/* write data fit into cache buffer, so return */
	}

	/* there is still data left to write - flush current buffer and load next block into cache */

	if((status = isioFlush(ifp)) < 0) /* flush the full cache buffer */
		return(status);

	/* mark cache buffer as empty */
	ifp->ptr = ifp->base;
	ifp->count = 0;

	wrtSize = countToWrite - wrtSize;
	if(wrtSize > ISIO_MAX_BUF_SIZE) {
		/* Too big to cache, just send it */
		dataObjWriteOutBBuf.buf = srcPtr;
		dataObjWriteOutBBuf.len = wrtSize;

		memset(&dataObjWriteInp, 0, sizeof(dataObjWriteInp));

		dataObjWriteInp.l1descInx = ifp->l1descInx;
		dataObjWriteInp.len = wrtSize;

		status = rcDataObjWrite(Comm, &dataObjWriteInp, &dataObjWriteOutBBuf);

		if(debug)
			printf("isioFileWrite: rcDataWrite 2 %d\n", status);

		if(status < 0)
			return(status);

		/* update base_offset */
		ifp->base_offset += (ifp->bufferSize + wrtSize);

		return(status);  /* total bytes written */
	}

	/* write the remaining user data to cache buffer */

	/* expand the cache buffer if necessary */
	if(wrtSize > ifp->bufferSize) {
		/* reallocate larger buffer */
		newBufSize=(2*wrtSize)+8;  /* Possible next size */
		if((ifp->base=realloc(ifp->base, newBufSize)) == NULL) {
			fprintf(stderr,"Memory Allocation error\n");
			return(0);
		}
		ifp->ptr = ifp->base;
	}

	memcpy(ifp->ptr, srcPtr, wrtSize);
	ifp->dirty = 1;
	ifp->ptr += wrtSize;
	ifp->count -= wrtSize;
	if(ifp->count < 0)
		ifp->count = 0;
	
	write_count += wrtSize;

	return(write_count);
}

size_t irodsfwrite(void *buffer, size_t itemsize, size_t nitems, FILE *fi_stream)
{
	int i;
	ISIO_FILE *ifp;

	i = (int)(unsigned long)fi_stream;
	if (debug) printf("irodsfwrite: %d\n", i);

	if((ifp=is_in_map(i)) != NULL) {
		return(isioFileWrite(ifp, buffer, itemsize*nitems)/itemsize);
	}
	else {
		return(fwrite(buffer, itemsize, nitems, fi_stream));
	}
}

int isioFileClose(ISIO_FILE *ifp)
{
	openedDataObjInp_t dataObjCloseInp;
	int status;

	if(debug) printf("isioFileClose\n");

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	/* If the buffer had been used for writing, flush it */
	if((status = isioFlush(ifp)) < 0) return(status);

	memset (&dataObjCloseInp, 0, sizeof (dataObjCloseInp));
	dataObjCloseInp.l1descInx = ifp->l1descInx;

	destroy_isio_file(ifp);

	return(rcDataObjClose(Comm, &dataObjCloseInp));
}

int irodsfclose(FILE *fi_stream)
{
	int i;
	ISIO_FILE *ifp;

	i = (int)(unsigned long)fi_stream;
	if (debug) printf("isiofclose: %d\n", i);

	if((ifp=is_in_map(i)) != NULL) {
		return(isioFileClose(ifp));
	}
	else {
		return(fclose(fi_stream));
	}
}

int isioFileTell(ISIO_FILE *ifp)
{
	openedDataObjInp_t seekParam;
	fileLseekOut_t* seekResult = NULL;
	int status;

	if (debug) printf("isioFileTell\n");

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	memset(&seekParam,  0, sizeof(openedDataObjInp_t));
	seekParam.l1descInx = ifp->l1descInx;
	seekParam.offset  = 0L;
	seekParam.whence  = SEEK_CUR;
	if((status = rcDataObjLseek(Comm, &seekParam, &seekResult )) < 0) {
		rodsLogError(LOG_ERROR, status, "isioFileSeek");
	}

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	return((uint64_t)(seekResult->offset));
}

long irodsftell(FILE *fi_stream)
{
	int i;
	ISIO_FILE *ifp;

	i = (int)(unsigned long)fi_stream;

	if (debug) printf("irodstell: %d\n", i);

	if((ifp=is_in_map(i)) != NULL) {
		return(isioFileTell(ifp));
	}
	else {
		return(ftell(fi_stream));
	}
}

/*
If the new read position (offset) is within the cached data, simply adjust cacheinfo values.
Otherwise, perform the seek using rcDataObjLseek() and set cacheinfo values to reflect absence of cached data
*/
int isioFileSeek(ISIO_FILE *ifp, long offset, int whence)
{
	openedDataObjInp_t seekParam;
	fileLseekOut_t* seekResult = NULL;
	int status = 0;
	long adj = 0L;

	if (debug) printf("isioFileSeek\n");

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	if(offset_in_current_buffer(ifp, offset, whence, &adj)) {
		ifp->ptr += adj;
		ifp->count -= adj;
	}
	else {
		if((status = isioFlush(ifp)) < 0)
			return(status);
		memset(&seekParam,  0, sizeof(openedDataObjInp_t));
		seekParam.l1descInx = ifp->l1descInx;
		seekParam.offset  = offset;
		seekParam.whence  = whence;
		if((status = rcDataObjLseek(Comm, &seekParam, &seekResult)) < 0) {
			rodsLogError(LOG_ERROR, status, "isioFileSeek");
		}

		ifp->base_offset = isioFileTell(ifp) - status;

		/* mark cache buffer as empty */
		ifp->ptr = ifp->base;
		ifp->count = 0;
	}

	if(status > 0)
		status = 0;	/* fseek should return 0 for success, so adjust */

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	return(status);
}

/*
offset_in_current_buffer:
determine if requested position is within the current cacheInfo buffer
	if it isn't, return 0
	otherwise: return 1 and return calculated adjustment in *adj (if adj is not NULL)
*/
static int offset_in_current_buffer(ISIO_FILE *ifp, long offset, int whence, long *adj)
{
	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	if(ifp->base == NULL || ifp->ptr == NULL || ifp->bufferSize <= 0 || (ifp->ptr - ifp->base) + ifp->count <= 0) {
		return(0);
	}

	if(whence == SEEK_SET) {

		if(offset >= 0
			&& offset >= ifp->base_offset
			&& offset < (ifp->base_offset + ifp->ptr - ifp->base) + ifp->count) {

			if(adj != NULL) {
				*adj = (offset - ifp->base_offset) - (ifp->ptr - ifp->base);
			}

			return(1);
		}
		else {
			return(0);
		}

	}
	else if(whence == SEEK_CUR) {

		if(offset > 0) {
			if(offset <= ifp->count) {
				if(adj != NULL)
					*adj = offset;

				return(1);
			}
			else {
				return(0);
			}
		}
		else if(offset < 0) {
			if(((ifp->ptr - ifp->base) + offset) > 0) {
				if(adj != NULL)
					*adj = offset;

				return(1);
			}
			else {
				return(0); 
			}
		}
		else {
			if(adj != NULL)
				*adj = 0;
			return(1);	/* zero offset must be in buffer */
		}

	}
	else if(whence == SEEK_END) {
		return(0);	/* just do the data read to buffer */
	}
	else {
		return(0);	/* unrecognised whence value ignored? */
	}

	return(0);
}

int irodsfseek(FILE *fi_stream, long offset, int whence)
{
	int i;
	ISIO_FILE *ifp;

	i = (int)(unsigned long)fi_stream;

	if (debug) printf("isiofseek: %d\n", i);

	if((ifp=is_in_map(i)) != NULL) {
		return(isioFileSeek(ifp, offset, whence));
	}
	else {
		return(fseek(fi_stream, offset, whence));
	}
}

/*
isioFlush:
	This function should ensure that the contents of the iRODS store match
	what's in the cache buffer. It doesn't handle adjustment of any of the
	cache buffer pointers or counts.
	if dirty flag is set, write the contents of the current cache buffer
		and clear dirty flag.
	Caller is responsible for any buffer freeing or (re)allocation.
	Returns: < 0 indicates error, >= 0 success
*/
int isioFlush(ISIO_FILE *ifp)
{
	openedDataObjInp_t dataObjWriteInp;
	bytesBuf_t dataObjWriteOutBBuf;
	int status;

	if (debug) printf("isioFlush\n");

	assert((ifp->base != NULL) && (ifp->ptr >= ifp->base) && ((ifp->ptr - ifp->base) + ifp->count <= ifp->bufferSize));

	if(ifp->dirty) {

		dataObjWriteOutBBuf.buf = ifp->base;
		dataObjWriteOutBBuf.len = (ifp->ptr - ifp->base) + ifp->count;

		memset(&dataObjWriteInp, 0, sizeof(dataObjWriteInp));

		dataObjWriteInp.l1descInx = ifp->l1descInx;
		dataObjWriteInp.len = dataObjWriteOutBBuf.len;	/* write out entire buffer instead of trying to track all writes to buffer */

		if(debug) printf("isioFlush: writing %d\n", dataObjWriteOutBBuf.len);
		status = rcDataObjWrite(Comm, &dataObjWriteInp, &dataObjWriteOutBBuf);
		if(status >= 0) {
			ifp->dirty = 0;
		}

		return(status);
	}

	return(0);
}

int irodsfflush(FILE *fi_stream)
{
	int i;
	ISIO_FILE *ifp;

	i = (int)(unsigned long)fi_stream;

	if (debug) printf("isiofflush: %d\n", i);

	if((ifp=is_in_map(i)) != NULL) {
		return(isioFlush(ifp));
	}
	else {
		return(fflush(fi_stream));
	}
}

int isioFilePutc(int inchar, ISIO_FILE *ifp)
{
	int mychar;

	mychar = inchar;

	return(isioFileWrite(ifp, &mychar, 1));
}

int irodsfputc(int inchar, FILE *fi_stream)
{
	int i;
	ISIO_FILE *ifp;

	i = (int)(unsigned long)fi_stream;
	if (debug) printf("isiofputc: %d\n", i);

	if((ifp=is_in_map(i)) != NULL) {
		return(isioFilePutc(inchar, ifp));
	}
	else {
		return(fputc(inchar, fi_stream));
	}
}

int isioFileGetc(ISIO_FILE *ifp)
{
	int mychar=0;
	int status;

	status = isioFileRead(ifp, &mychar, 1);

	if (status==0) return(EOF);

	if (status<0) return(status);

	return(mychar);
}

int irodsfgetc(FILE *fi_stream)
{
	int i;
	ISIO_FILE *ifp;

	i = (int)(unsigned long)fi_stream;
	if (debug) printf("isiofgetc: %d\n", i);

	if((ifp=is_in_map(i)) != NULL) {
		return(isioFileGetc(ifp));
	}
	else {
		return(fgetc(fi_stream));
	}
}


void irodsexit(int exitValue)
{
	int status;

	if (debug) printf("irodsexit: %d\n", exitValue);

	if (setupFlag>0) {
		status = rcDisconnect(Comm);
	}
	exit(exitValue);
}

