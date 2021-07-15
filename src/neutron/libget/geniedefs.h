#ifndef ISISLIBGET
#define ISISLIBGET
/*
 * $Id: libget.h,v 1.5 1999/02/04 13:36:54 faa Exp $
 *
 * Uncomment the definition of WORDS_BIGENDIAN if you are on a BIG ENDIAN machine
 * i.e. a machine where the most significant byte is stored first in memory
 *
 * Examples of BIG ENDIAN machines are SUN, HP, Silicon Graphics and Apple Macintosh
 */

/* #define WORDS_BIGENDIAN 	1 */


/*
 * Uncomment this definition if you are using f2c or g77, or get a lot
 * of undefined symbols on making the example program
 */

/* #define HAVE_LIBF2C		1 */


/*
 * Uncomment this only if you have the DCE (Distributed Computing Environment)
 * software installed and you would also like to read directly from an ISIS DAE
 */

/* #define HAVE_LIBDCE		1 */

/* Uncomment this if you are using GCC under MinGW. */

#define GCC_MINGW          1 

/* Uncomment this if you are using g95 with GCC under MinGW */

/* #define G95_MINGW          1 */

/* 
 * leave alone after this point 
 */

#if (defined(__unix__) || defined(__APPLE__)) && !defined(__unix)
#define __unix 1
#endif /* defined(__unix__) && !defined(__unix) */

#ifdef _WIN32
   #if (defined(GCC_MINGW) || defined(G95_MINGW))
      #define GENIEMECHANISM __cdecl
   #else
      #define GENIEMECHANISM __stdcall
   #endif
#else
   #define GENIEMECHANISM
#endif  


#define USE_RCSID(rcsid)        static void* use_##rcsid = (void*)&rcsid
#define g_printef gprintef

extern void g_printef(const char *fmt, ...);
extern void gprintif(const char *fmt, ...);
extern void gprintdf(const char *fmt, ...);

#endif /* ISISLIBGET */
