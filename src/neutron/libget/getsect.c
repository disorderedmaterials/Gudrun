#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/*
 * File: 	getsect.c - A fast interface to the old GenieII get routines.
 * Author:	Freddie Akeroyd, ISIS computer group
 *
 * Basically, we just read the whole file into "file_buffer" and then catch
 * calls to the old FORTRAN GETSECT routine
 */
#include "geniedefs.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <errno.h>
#if HAVE_FCNTL_H
#include <fcntl.h>
#endif /* HAVE_FCNTL_H */

static const char* rcsid = "@(#)$Id: getsect.c,v 1.18 1999/02/04 11:51:16 faa Exp $";
USE_RCSID(rcsid);

#ifdef HAVE_UNISTD_H
#	include <unistd.h>
#endif

#define AD_TCB	26	/* from crptsect.def */

#ifdef HAVE_LIBDCE
#include <dae/rio_interface.h>			/* DAE reading routines */
#else
#define RIO_OPERATION_FAILED 			0
#define RIO_OPERATION_SUCCESSFUL        	1
#define CLOSE_DAE()				RIO_OPERATION_FAILED
#define open_dae_cds(__file)			RIO_OPERATION_SUCCESSFUL	/* hack for libget */
#define RMEMVI(__a, __b, __c, __d, __e) 	RIO_OPERATION_FAILED
#endif /* HAVE_LIBDCE */

#ifdef __VMS_POSIX
#	include <vms_posix_extensions.h>	/* for _translate() */
#endif
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
#	include <io.h>
#	include <fcntl.h>
#else
#   include <unistd.h>
#	include <sys/file.h>
#endif /* _WIN32 */


#if !defined(USE_OPEN) && !defined(USE_FOPEN)
#define USE_FOPEN
#endif

/* Sort out external symbol names */
#if defined(HAVE_LIBF2C)
#define fastget_init fastget_init__
#define getsect getsect_
#define fort_file fort_file__
#elif defined(__unix)
#define fastget_init fastget_init_
#define getsect getsect_
#define fort_file fort_file_
#elif defined(_WIN32)
   #if defined(GCC_MINGW)
      #define fastget_init fastget_init_
      #define getsect getsect_
      #define fort_file fort_file_
   #elif defined(G95_MINGW)
      #define fastget_init fastget_init__
      #define getsect getsect_
      #define fort_file fort_file__
   #else
      #define fastget_init FASTGET_INIT
      #define getsect GETSECT
      #define fort_file FORT_FILE
   #endif
#endif

int GENIEMECHANISM fastget_init(const char* name, const int* iunit);
void GENIEMECHANISM getsect(const int* istart, const int* ilong, int ivalue[], const int* iunit, int* ierr);
void GENIEMECHANISM fort_file(const int* iunit, char* file_name);
void error_add(const char *obj, const char *fail, const char *soln);

static int dae_access = 0;					/* are we reading from ISIS DAE? */
static int crpt_access = 0;					/* are we reading from ISIS CRPT? */
static int start_of_data_section = 0;				/* offset of data section (words) in CRPT i.e. the */
								/* point after which we want to query the DAE */
static char* crpt_address = NULL;				/* pointer to global section */
static long crpt_size = 0;					/* crpt size in words */
static int dae_word_length = 0;

static int enlarge_buffer(off_t new_buffer_size);

typedef enum { ReadWholeFile, UseFOPEN } FastgetBufferOption;

static char file_name[257];					/* Name of the file to buffer (passed from fortran) */
static char* file_buffer = NULL;				/* Where we store "file_name" */
static off_t file_size = 0;					/* How big "file_name" is */
static off_t file_buffer_size = 0;				/* How big our buffer is ( file_size <= file_buffer_size ) */
static long io_block_size = 0;					/* Optimal read block size for the file system containing file_name */
static long buffer_iunit = 0;					/* Fortran unit associated with "file_buffer" */
static FastgetBufferOption buffer_option = ReadWholeFile;	/* how "file buffer" is being used */
/* 
 * G_MAX_BUFFER_SIZE is the maximum size the buffer may be
 * G_FOPEN_BUFFER_SIZE is the buffer size used if file is larger than G_MAX_BUFFER_SIZE
 */
static off_t MAX_BUFFER_SIZE = 0;
static off_t FOPEN_BUFFER_SIZE = 0;
static FILE* file_fd = NULL;

/* 
 * Returns 1 on error, 0 on success
 * name is character*120 fortran variable 
 */

#ifdef __VMS
#include <ssdef.h>
#include <secdef.h>
#include <stsdef.h>
#include <lib$routines.h>
#include <descrip.h>
#include <starlet.h>
#include <unixio.h>
#include <fcntl.h>

/* Modify from Severe error to error, otherwise program will exit in call to lib$signal() */
#define VMS_FATAL_TO_ERROR(__status) \
        if ( (__status & 07) == STS$K_SEVERE ) \
        { \
            __status &= ~STS$K_SEVERE;  \
            __status |= STS$K_ERROR;  \
        }


/*
 *  This routines maps to global section "name", returning the memory 
 *  addresses in ret_addr[].
 *  return 0 on success, 1 on error
 */
static int map_global(const char *name, int section, int* ret_addr[])
{
    int* in_addr[2];
    unsigned status, flags=SEC$M_SYSGBL|SEC$M_EXPREG;
    struct dsc$descriptor_s section_name;
    in_addr[0] = in_addr[1] = (int*)200;  /* force extra space into P0 */

    section_name.dsc$w_length = strlen(name);
    section_name.dsc$a_pointer = (char*)name;
    section_name.dsc$b_class = DSC$K_CLASS_S;
    section_name.dsc$b_dtype = DSC$K_DTYPE_T;

    status = sys$mgblsc(in_addr, ret_addr, 0, flags, &section_name, 0, 0);
    if (status != SS$_NORMAL)
    {
	gprintef("MAP_GLOBAL: Error code %u returned from sys$mgblsc() while mapping to \"%s\"\n", status, name);
	VMS_FATAL_TO_ERROR(status);
	lib$signal(status);
	return 1;
    }
    return 0;
}

#else
static int map_global(const char *name, int section, int* ret_addr[])
{
#ifdef HAVE_LIBDCE
    int section_size;
    if (READ_CRPT(&ret_addr[0], &section, &section_size) == RIO_OPERATION_SUCCESSFUL)
    {
        ret_addr[1] = ret_addr[0] + section_size;
        return 0;
    }
    else
    {
	ret_addr[1] = ret_addr[0] = 0;
        gprintef("GETSECT: DAE Access not available\n");
        return 1;
    }
#else
    gprintef("GETSECT: DAE Access not available\n");
    return 1;
#endif /* HAVE_LIBDCE */
}
#endif /* __VMS */

/* 0 on success, 1 on failure */
int GENIEMECHANISM fastget_init(const char* name, const int* iunit)
{
    static char error_message[256], cds_name[256];
    struct stat stat_buffer;
    static int first_call = 1;
/*    stat_t stat_buffer; */
    int i, ilen, *icrpt;
    int daepoint;
    int* ret_addr[2];
    char sect_name[32];
    const char* inst_name = getenv("INST_NAME");
    char* temp_str;
    if (inst_name == NULL)
    {
	inst_name = "unknown";
    }
    if (first_call == 1)
    {
	first_call = 0;
	MAX_BUFFER_SIZE = (getenv("G_MAX_BUFFER_SIZE") != 0 ? atol(getenv("G_MAX_BUFFER_SIZE")) : 4000000);
	FOPEN_BUFFER_SIZE = (getenv("G_FOPEN_BUFFER_SIZE") != 0 ? atol(getenv("G_FOPEN_BUFFER_SIZE")) : BUFSIZ);
    }
/* First null terminate our FORTRAN string */
/* Modified 12/22 by TGAY - Change loop limit from 119 to 255 for consistency with string variable declared size */
    for(i=255; name[i] == ' '; i--)
	;
    if (i < 0)
    {
	error_add("FASTGET_INIT", "null file name given", "");
	return 1;
    }
    strncpy(file_name, name, i+1);
    file_name[i+1] = '\0';
#ifdef __VMS_POSIX
    const char* trans_file_name = _translate(file_name, _TO_PATHNAME);
    if (trans_file_name != NULL)
    {
	strcpy(file_name, trans_file_name);
    }
#endif /* __VMS_POSIX */
    buffer_iunit = *iunit;
/* Now attempt to open the file and get some information about it */
    dae_access = crpt_access = crpt_size = 0;
    crpt_address = NULL;
    if ( (strncmp(file_name, "/.../", 5) == 0) || (strncmp(file_name, "/.:/", 4) == 0) )
    {
	i = strlen(file_name);
	strcpy(cds_name, file_name);
	if ( strcmp(file_name + i - 4, "_dae") == 0 )
	{
	    dae_access = 1;
	}
	else if ( strcmp(file_name + i - 5, "_crpt") == 0 )
	{
	    crpt_access = 1;
	    strcpy(cds_name + i - 5, "_dae");	/* we only export _dae to CDS, so fudge name */
	}
    }
    if (dae_access || crpt_access)
    {
/* 
 * name is of the form /.:/servers/inst_dae or /.:/servers/inst_crpt
 */
	temp_str = strstr(cds_name, "/servers/");
	if (temp_str == NULL)
	{
	    sprintf(error_message, "DAE name \"%s\" does not follow a recognised format", cds_name); 
	    error_add("FASTGET_INIT", error_message, "");
	    return 1;
	}
	CLOSE_DAE();
	if (open_dae_cds(cds_name) == RIO_OPERATION_FAILED)
	{
	    sprintf(error_message, "Failed to open DAE at \"%s\"", cds_name); 
	    error_add("FASTGET_INIT", error_message, "");
	    return 1;
	}
	temp_str += strlen("/servers/");	/* move past /servers/ */
	ilen = strlen(temp_str) - 4;	/* to knock off "_dae" suffix */
	strcpy(sect_name, "G_");
	for(i=0; i<ilen; i++)
	{
	    sect_name[2+i] = toupper(temp_str[i]);
	}
	sect_name[ilen+2] = '\0';
	i = (dae_access ? AD_TCB : 0);
	if (map_global(sect_name, i, ret_addr) != 0)
	{
	    return 1;
	}
	crpt_address = (char*)ret_addr[0];
	crpt_size = ret_addr[1] - ret_addr[0] + 1;
	gprintdf("GETSECT: Mapping to portion of CRPT size = %d words\n", crpt_size);
	icrpt = ret_addr[0];
	if (icrpt[20] == 1)	/* test VER1 */
	{
	    start_of_data_section = icrpt[26];	/* DATAPOINT - section addresses start at offset 21 */
	}
	else
	{
	    start_of_data_section = icrpt[27];	/* DATAPOINT - section addresses start at offset 21 */
	}
	daepoint = icrpt[24];	/* DAE section of CRPT, IFORMAT(4) in io.f */
	if (daepoint < crpt_size)
	{
	    dae_word_length = icrpt[daepoint];
	    gprintdf("GETSECT: DAE word length = %d bytes\n", dae_word_length);
	}
	else
	{
	    gprintef("GETSECT: error in getting DAE word length\n");    
	}
	gprintdf("GETSECT: data section starts at %d words\n", start_of_data_section);
    }
    else
    {
#if defined(__VMS) && !defined(__VMS_POSIX)
        int fd = open(file_name, O_RDONLY, 0, "shr=put");
#else
        int fd = open(file_name, O_RDONLY, 0);
#endif /* __VMS && !__VMS_POSIX */
        if (fd < 0) 
        {
	    sprintf(error_message, "Cannot open \"%s\" --- %s", file_name, strerror(errno));
	    error_add("FASTGET_INIT", error_message, "Check file exists and is readable");
	    return 1;
        }
/***** ???   if (stat(file_name, &stat_buffer) != 0)   seemed to fail sometimes ??? *****/
        if (fstat(fd, &stat_buffer) != 0)
        {
	    sprintf(error_message, "Cannot fstat \"%s\" --- %s", file_name, strerror(errno));
	    error_add("FASTGET_INIT", error_message, "Check file is not corrupted");
            close(fd);
	    return 1;
        }
        if (close(fd) != 0)
        {
	    sprintf(error_message, "Cannot close \"%s\" --- %s", file_name, strerror(errno));
	    error_add("FASTGET_INIT", error_message, "This shouldn't happen!");
	    return 1;
        }
        file_size = stat_buffer.st_size;
#if defined(__VMS) || defined(_WIN32)
        io_block_size = 512;			/* optimal IO block size for this file system */
#else
        io_block_size = stat_buffer.st_blksize;	/* optimal IO block size for this file system */
#endif /* __VMS */
        if (file_size > MAX_BUFFER_SIZE)
        {
	    buffer_option = UseFOPEN;
	    enlarge_buffer(FOPEN_BUFFER_SIZE);
        }
        else
        {
	    buffer_option = ReadWholeFile;
	    if (file_size > file_buffer_size)
	    {
	        enlarge_buffer(file_size);
	    }
        }

#ifdef USE_FOPEN
        if (file_fd != 0)
        {
	    fclose(file_fd);
	    file_fd = 0;
        }
/* Now fill up our buffer */
#if defined(__VMS) && !defined(__VMS_POSIX)
        if ((file_fd = fopen(file_name, "r", "shr=put", "ctx=stm")) == NULL)
#else
        if ((file_fd = fopen(file_name, "rb")) == NULL)
#endif /* __VMS && !__VMS_POSIX */
        {
	    sprintf(error_message, "Cannot fopen \"%s\" --- %s", file_name, strerror(errno));
	    error_add("FASTGET_INIT", error_message, "Check file exists and is readable");
	    return 1;
        }
        if ((buffer_option == ReadWholeFile) && (fread(file_buffer, 1, file_size, file_fd) != (size_t)file_size))
        {
	    sprintf(error_message, "Cannot read %ld bytes from \"%s\" --- %s", (long)file_size, file_name, strerror(errno));
	    error_add("FASTGET_INIT", error_message, "Check file is not corrupted or has been deleted");
	    return 1;
        }          
        else if (buffer_option == UseFOPEN)
        {
	    if (setvbuf(file_fd, file_buffer, _IOFBF, FOPEN_BUFFER_SIZE) != 0)
	    {
	        perror("setvbuf");
	        return 1;
	    }
        }
#endif /* USE_FOPEN */
#ifdef USE_STREAMS
        ifstream the_file;
        the_file.open(file_name);
        the_file.read(file_buffer, file_size);
        the_file.close();
#endif /* USE_STREAMS */
    }
    return 0;
}

/*
 * Originally declared as SUBROUTINE GETSECT(ISTART,ILONG,IVALUE,IUNIT,IERR) in FORTRAN
 * 
 * The routine read a certain number of 512 byte blocks from a direct access file on unit "iunit"
 */

static void getsect_crpt(const int* istart, const int* ilong, int ivalue[], const int* iunit, int* ierr)
{
    static const int sz = sizeof(int);
    static char buffer[256];
    short zbdev = 0;
    *ierr = 0;
    if ( (*ilong > 0) && ((*istart + *ilong - 1) <= crpt_size) )
    {
    	memcpy(ivalue, crpt_address + (*istart - 1)*sz, (*ilong) * sz);
    }
    else
    {
	sprintf(buffer, "Attempt to read invalid part of CRPT buffer: istart = %d words, ilong = %d words, CRPT size = %ld words",
			*istart, *ilong, crpt_size);
	error_add("GETSECT", buffer, "This shouldn't happen!");
	*ierr = -1;
	return;
    }
}

static void getsect_dae(const int* istart, const int* ilong, int ivalue[], const int* iunit, int* ierr)
{
    int start_of_spectra;
    static const int sz = sizeof(int);
    static char buffer[256];
    short zbdev = 0;
    *ierr = 0;
/* We want to read everything but the counts from the CRPT, so check the range */
    if ( *istart > start_of_data_section )
    {
/* need to DAE read here */
	start_of_spectra = (*istart - start_of_data_section - 1) * dae_word_length;
	if (RMEMVI(&start_of_spectra, ivalue, (fort_i4*)ilong, &zbdev, &dae_word_length) == RIO_OPERATION_FAILED)
	{
	    sprintf(buffer, "Attempt to read %d words from DAE (word length %d) staring at address %d failed",
			*ilong, dae_word_length, start_of_spectra);
	    error_add("GETSECT", buffer, "Is DAE access supported on this platform?");
	}
    }
    else if ( (*ilong > 0) && ((*istart + *ilong - 1) <= crpt_size) )
    {
    	memcpy(ivalue, crpt_address + (*istart - 1)*sz, (*ilong) * sz);
    }
    else
    {
	sprintf(buffer, "Attempt to read invalid part of CRPT buffer: istart = %d words, ilong = %d words, CRPT size = %ld words",
			*istart, *ilong, crpt_size);
	error_add("GETSECT", buffer, "This shouldn't happen!");
	*ierr = -1;
	return;
    }
}

static void getsect_file(const int* istart, const int* ilong, int ivalue[], const int* iunit, int* ierr)
{
    static char buffer[256];
    static const int sz = sizeof(int);
#if 0
    char fort_file_name[120];		/* fortran file attached to unit "iunit" */
    if (*iunit != buffer_iunit)		/* If we've changed file, reinitialize */
    {
	fort_file(*iunit, fort_file_name);
	fastget_init(fort_file_name, *iunit);
    }
#endif
    if ((*ilong < 0) || ((*istart + *ilong - 1) * sz > file_size))
    {
	sprintf(buffer, "Attempt to read invalid part of file buffer: istart = %d words, ilong = %d words, file size = %ld bytes",
			*istart, *ilong, file_size);
	error_add("GETSECT", buffer, "This shouldn't happen!");
	*ierr = -1;
	return;
    }
    *ierr = 0;
    if (buffer_option == ReadWholeFile)
    {
    	memcpy(ivalue, file_buffer + (*istart - 1)*sz, (*ilong) * sz);
    }
    else
    {
	if ( fseek(file_fd, (*istart - 1)*sz, SEEK_SET) == -1 )
	{
	    sprintf(buffer, "Error seeking file size %ld bytes at %d: %s", 
			file_size, (*istart - 1)*sz, strerror(errno));
	    error_add("GETSECT", buffer, "This shouldn't happen!");
    	    *ierr = -1;
	    return;
	}
        if ( fread(ivalue, sz, *ilong, file_fd) != (size_t)(*ilong) )
        {
	    sprintf(buffer, "Error reading file size %ld bytes for %d words from %d: %s", 
			file_size, *ilong, *istart - 1, strerror(errno));
	    error_add("GETSECT", buffer, "This shouldn't happen!");
	    *ierr = -1;
	    return;
	}
    }
    return; 
}

void GENIEMECHANISM getsect(const int* istart, const int* ilong, int ivalue[], const int* iunit, int* ierr)
{
    if (dae_access)
    {
	getsect_dae(istart, ilong, ivalue, iunit, ierr);
    }
    else if (crpt_access)
    {
	getsect_crpt(istart, ilong, ivalue, iunit, ierr);
    }
    else
    {
	getsect_file(istart, ilong, ivalue, iunit, ierr);
    }
    return; 
}

static int enlarge_buffer(off_t new_buffer_size)
{
    static char error_message[256];
    if (new_buffer_size > file_buffer_size)
    {
	if (file_buffer != NULL)
	    free(file_buffer);
	file_buffer = (char*)malloc(new_buffer_size * sizeof(char));
	if (file_buffer == NULL)
	{
	    sprintf(error_message, "Unable to allocate %ld bytes for file buffer", (long)file_size);
	    error_add("FASTGET_INIT", error_message, "Check your page file quota and/or the current machine load");
	    file_buffer_size = 0;
	    return 1;
	}
	file_buffer_size = new_buffer_size;
    }        
    return 0;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

