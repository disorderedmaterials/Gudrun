#ifndef INT_CONVERT
#define INT_CONVERT

/* 
 * $Id: int_convert.h,v 1.7 1999/01/07 13:10:39 faa Exp $
 */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* to allow calling from fortran */
#if defined(HAVE_LIBF2C)
#define local_to_vax_short local_to_vax_short__
#define vax_to_local_short vax_to_local_short__
#define local_to_vax_int local_to_vax_int__
#define vax_to_local_int vax_to_local_int__
#define local_to_vax_shorts local_to_vax_shorts__
#define vax_to_local_shorts vax_to_local_shorts__
#define local_to_vax_ints local_to_vax_ints__
#define vax_to_local_ints vax_to_local_ints__
#elif defined(__unix)
#define local_to_vax_short local_to_vax_short_
#define vax_to_local_short vax_to_local_short_
#define local_to_vax_int local_to_vax_int_
#define vax_to_local_int vax_to_local_int_
#define local_to_vax_shorts local_to_vax_shorts_
#define vax_to_local_shorts vax_to_local_shorts_
#define local_to_vax_ints local_to_vax_ints_
#define vax_to_local_ints vax_to_local_ints_
#elif defined(_WIN32)
   #if defined(GCC_MINGW) 
      #define local_to_vax_short local_to_vax_short_
      #define vax_to_local_short vax_to_local_short_
      #define local_to_vax_int local_to_vax_int_
      #define vax_to_local_int vax_to_local_int_
      #define local_to_vax_shorts local_to_vax_shorts_
      #define vax_to_local_shorts vax_to_local_shorts_
      #define local_to_vax_ints local_to_vax_ints_
      #define vax_to_local_ints vax_to_local_ints_
   #elif defined(G95_MINGW)
      #define local_to_vax_short local_to_vax_short__
      #define vax_to_local_short vax_to_local_short__
      #define local_to_vax_int local_to_vax_int__
      #define vax_to_local_int vax_to_local_int__
      #define local_to_vax_shorts local_to_vax_shorts__
      #define vax_to_local_shorts vax_to_local_shorts__
      #define local_to_vax_ints local_to_vax_ints__
      #define vax_to_local_ints vax_to_local_ints__
   #else
      #define local_to_vax_short LOCAL_TO_VAX_SHORT
      #define vax_to_local_short VAX_TO_LOCAL_SHORT
      #define local_to_vax_int LOCAL_TO_VAX_INT
      #define vax_to_local_int VAX_TO_LOCAL_INT
      #define local_to_vax_shorts LOCAL_TO_VAX_SHORTS
      #define vax_to_local_shorts VAX_TO_LOCAL_SHORTS
      #define local_to_vax_ints LOCAL_TO_VAX_INTS
      #define vax_to_local_ints VAX_TO_LOCAL_INTS
   #endif
#endif 


unsigned short GENIEMECHANISM local_to_vax_short(const unsigned short* s);
unsigned short GENIEMECHANISM vax_to_local_short(const unsigned short* s);
unsigned GENIEMECHANISM local_to_vax_int(const unsigned* i);
unsigned GENIEMECHANISM vax_to_local_int(const unsigned* i);

void GENIEMECHANISM local_to_vax_shorts(unsigned short* sa, const int* n);
void GENIEMECHANISM vax_to_local_shorts(unsigned short* sa, const int* n);
void GENIEMECHANISM local_to_vax_ints(unsigned* ia, const int* n);
void GENIEMECHANISM vax_to_local_ints(unsigned* ia, const int* n);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* INT_CONVERT */
