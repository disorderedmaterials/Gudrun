#ifndef FLOAT_CONVERT
#define FLOAT_CONVERT

/* 
 * $Id: float_convert.h,v 1.7 1999/01/07 13:10:38 faa Exp $ 
 */

/* this is to allow the routines to be called from fortran */
#if defined(HAVE_LIBF2C)
#define ieee_to_vax_float ieee_to_vax_float__
#define vax_to_ieee_float vax_to_ieee_float__
#define vaxf_to_local vaxf_to_local__
#define local_to_vaxf local_to_vaxf__
#elif defined(__unix)
#define ieee_to_vax_float ieee_to_vax_float_
#define vax_to_ieee_float vax_to_ieee_float_
#define vaxf_to_local vaxf_to_local_
#define local_to_vaxf local_to_vaxf_
#elif defined(_WIN32)
   #if defined(GCC_MINGW)
      #define ieee_to_vax_float ieee_to_vax_float_
      #define vax_to_ieee_float vax_to_ieee_float_
      #define vaxf_to_local vaxf_to_local_
      #define local_to_vaxf local_to_vaxf_
   #elif defined(G95_MINGW)
      #define ieee_to_vax_float ieee_to_vax_float__
      #define vax_to_ieee_float vax_to_ieee_float__
      #define vaxf_to_local vaxf_to_local__
      #define local_to_vaxf local_to_vaxf__
   #else
      #define ieee_to_vax_float IEEE_TO_VAX_FLOAT
      #define vax_to_ieee_float VAX_TO_IEEE_FLOAT
      #define vaxf_to_local VAXF_TO_LOCAL
      #define local_to_vaxf LOCAL_TO_VAXF
   #endif
#endif

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* these routines return 0 = success, 1 = failure */

/* convert an IEEE single float to VAX F FLOAT format */
/* int ieee_to_vax_float(float* fp); */

/* convert VAX F FLOAT into an IEEE single float */
/* int vax_to_ieee_float(float* fp); */

/* convert float array val[n] to and from vax float */
void GENIEMECHANISM vaxf_to_local(float *val, const int *n, int *errcode);
void GENIEMECHANISM local_to_vaxf(float *val, const int *n, int *errcode);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* FLOAT_CONVERT */
