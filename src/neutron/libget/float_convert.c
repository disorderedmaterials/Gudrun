/*
 *      Module:         float_convert
 *      Author:         Freddie Akeroyd, ISIS
 *      Purpose:        Routines to convert from vax to ieee floating point, based on the XDR routines of SUN RPC
 *
 */
#include "geniedefs.h"
#include <stdio.h>
#include "float_convert.h"
#include "int_convert.h"

static const char* rcsid = "@(#)$Id: float_convert.c,v 1.9 1998/02/03 18:10:31 cmm Exp $";
USE_RCSID(rcsid);

/*
 * determine a few things we need to know to write machine independent data 
 */
 
#ifndef __VMS
#define IEEEFP			1
#endif /* __VMS */

#ifdef __VMS
/* set up codes for use of CVT$CONVERT_FLOAT */
#  if __IEEE_FP
#    define IEEEFP		1
#    define VMS_FLOAT_NATIVE	CVT$K_IEEE_S
#    define VMS_DOUBLE_NATIVE	CVT$K_IEEE_D
#  elif __D_FLOAT
#    define VAXFP		1
#    define VMS_FLOAT_NATIVE	CVT$K_VAX_F
#    define VMS_DOUBLE_NATIVE	CVT$K_VAX_D
#  elif __G_FLOAT
#    define VAXFP		1
#    define VMS_FLOAT_NATIVE	CVT$K_VAX_F
#    define VMS_DOUBLE_NATIVE	CVT$K_VAX_G
#  else
#    error Cannot determine VMS floating point format
#  endif
#endif /* __VMS */

#if WORDS_BIGENDIAN

/* What IEEE single precision floating point looks like on local machine */
struct	ieee_single {
	unsigned int	sign    : 1;
	unsigned int	exp     : 8;
	unsigned int	mantissa: 23;
};

/* Vax single precision floating point */
struct	vax_single {
	unsigned int	mantissa2 : 16;
	unsigned int	sign      : 1;
	unsigned int	exp       : 8;
	unsigned int	mantissa1 : 7;
};

#else

/* What IEEE single precision floating point looks like on local machine */
struct	ieee_single {
	unsigned int	mantissa: 23;
	unsigned int	exp     : 8;
	unsigned int	sign    : 1;
};

/* Vax single precision floating point */
struct	vax_single {
	unsigned int	mantissa1 : 7;
	unsigned int	exp       : 8;
	unsigned int	sign      : 1;
	unsigned int	mantissa2 : 16;
};

#endif /* WORDS_BIGENDIAN */

#define VAX_SNG_BIAS	0x81
#define IEEE_SNG_BIAS	0x7f

static const struct sgl_limits_struct {
	struct vax_single s;
	struct ieee_single ieee;
} sgl_limits[2] = {
	{{ 0x7f, 0xff, 0x0, 0xffff },	/* Max Vax */
	{ 0x0, 0xff, 0x0 }},		/* Max IEEE */
	{{ 0x0, 0x0, 0x0, 0x0 },	/* Min Vax */
	{ 0x0, 0x0, 0x0 }}		/* Min IEEE */
};

#define mmax sgl_limits[0]
#define mmin sgl_limits[1]

#if WORDS_BIGENDIAN

/* What IEEE double precision floating point looks like */
struct	ieee_double {
	unsigned int	mantissa2 : 32;
	unsigned int	sign      : 1;
	unsigned int	exp       : 11;
	unsigned int	mantissa1 : 20;
};

/* Vax double precision floating point */
struct  vax_double {
	unsigned int	mantissa4 : 16;
	unsigned int	mantissa3 : 16;
	unsigned int	mantissa2 : 16;
	unsigned int	sign      : 1;
	unsigned int	exp       : 8;
	unsigned int	mantissa1 : 7;
};

#else

/* What IEEE double precision floating point looks like */
struct	ieee_double {
	unsigned int	mantissa1 : 20;
	unsigned int	exp       : 11;
	unsigned int	sign      : 1;
	unsigned int	mantissa2 : 32;
};

/* Vax double precision floating point */
struct  vax_double {
	unsigned int	mantissa1 : 7;
	unsigned int	exp       : 8;
	unsigned int	sign      : 1;
	unsigned int	mantissa2 : 16;
	unsigned int	mantissa3 : 16;
	unsigned int	mantissa4 : 16;
};

#endif /* WORDS_BIGENDIAN */


#define VAX_DBL_BIAS	0x81
#define IEEE_DBL_BIAS	0x3ff
#define MASK(nbits)	((1 << nbits) - 1)

static struct dbl_limits {
	struct	vax_double d;
	struct	ieee_double ieee;
} dbl_limits[2] = {
	{{ 0x7f, 0xff, 0x0, 0xffff, 0xffff, 0xffff },	/* Max Vax */
	{ 0x0, 0x7ff, 0x0, 0x0 }},			/* Max IEEE */
	{{ 0x0, 0x0, 0x0, 0x0, 0x0, 0x0},		/* Min Vax */
	{ 0x0, 0x0, 0x0, 0x0 }}				/* Min IEEE */
};

/* VAX is little endian, so we may need to flip */
static int maybe_flip_bytes(void* p, size_t n)
{
#if WORDS_BIGENDIAN
    unsigned i;
    unsigned char c_tmp, *c_p = (unsigned char*)p;
    for(i=0; i<n/2; i++)
    {
	c_tmp = c_p[i];
	c_p[i] = c_p[n-i-1];
	c_p[n-i-1] = c_tmp;
    }
#endif /* WORDS_BIGENDIAN */
    return 0;
}
    
/* convert VAX F FLOAT into a local IEEE single float */
static int vax_to_ieee_float(float* fp)
{
    struct ieee_single is;
    struct vax_single vs;
    struct sgl_limits_struct;
    maybe_flip_bytes(fp, sizeof(float));
	vs = *((struct vax_single *)fp);
                switch(vs.exp){
                case 0 :
                        /* all vax float with zero exponent map to zero */
                        is = mmin.ieee ;
                        break ;
                case 2 :
                case 1 :
                        /* These will map to subnormals */
                        is.exp = 0 ;
                        is.mantissa = (vs.mantissa1 << 16) | vs.mantissa2;
                        /* lose some precision */
                        is.mantissa >>= 3 - vs.exp ;
                        is.mantissa += (1 << (20 + vs.exp)) ;
                        break ;
                case 0xff : /* mmax.s.exp */
                        if( vs.mantissa2 == mmax.s.mantissa2
                                && vs.mantissa1 == mmax.s.mantissa1)
                        {
                                /* map largest vax float to ieee infinity */
                                is = mmax.ieee ;
                                break ;
                        } /* else, fall thru */
                default :
                        is.exp = vs.exp - VAX_SNG_BIAS + IEEE_SNG_BIAS;
                        is.mantissa = (vs.mantissa1 << 16) | vs.mantissa2;
                }

                is.sign = vs.sign;
        *fp = *((float*)&is);
        return 0;
}

/* convert a local IEEE single float to little endian VAX F FLOAT format */
static int ieee_to_vax_float(float* fp)
{
    struct ieee_single is;
    struct vax_single vs;
    struct sgl_limits_struct;
    is = *((struct ieee_single*)fp);
                switch(is.exp) {
                case 0 :
                        if(is.mantissa == mmin.ieee.mantissa)
                        {
                                vs = mmin.s ;
                        } else {
                                unsigned tmp = is.mantissa >> 20 ;
                                if(tmp >= 4) {
                                        vs.exp = 2 ;
                                } else if (tmp >= 2) {
                                        vs.exp = 1 ;
                                } else {
                                        vs = mmin.s ;
                                        break ;
                                } /* else */
                                tmp = is.mantissa - (1 << (20 + vs.exp )) ;
                                tmp <<= 3 - vs.exp ;
                                vs.mantissa2 = tmp ;
                                vs.mantissa1 = (tmp >> 16);
                        }
                        break ;
                case 0xfe :
                case 0xff :
                        vs = mmax.s ;
                        break ;
                default :
                        vs.exp = is.exp - IEEE_SNG_BIAS + VAX_SNG_BIAS;
                        vs.mantissa2 = is.mantissa;
                        vs.mantissa1 = (is.mantissa >> 16);
                }

                vs.sign = is.sign;
    *fp = *((float*)&vs);
    maybe_flip_bytes(fp, sizeof(float));	/* Make little endian */
    return 0;
}

void GENIEMECHANISM vaxf_to_local(float *val, const int *n, int *errcode)
{
    int i;
    *errcode=0;
#if defined(VAXFP)
    /* nothing required */
#elif defined(IEEEFP)
        for(i=0; i<*n; i++)
        {
                if (vax_to_ieee_float(i+val) != 0)
                {
                        *errcode=1;
                }
        }
#else
#error Unknown floating point format
#endif
    return;
}

void GENIEMECHANISM local_to_vaxf(float *val, const int *n, int *errcode)
{
        int i;
        *errcode=0;
#if defined(VAXFP)
        /* nothing required */
#elif defined(IEEEFP)
        for(i=0; i<*n; i++)
        {
                if (ieee_to_vax_float(i+val) != 0)
                {
                        *errcode=1;
                }
        }
#else
#error Unknown floating point format
#endif
        return;
}

static unsigned char flip_bits(unsigned char cc)
{
    static int init = 0;
    static unsigned char ct[256];
    if (!init)
    {
        unsigned _i,_j;
        for(_i=0; _i<256; _i++)
        {
	    ct[_i] = 0;
 	    for(_j=0; _j<8; _j++)
	    {
    	        if (_i & (1<<_j)) { ct[_i] |= (128 >> _j); }
 	    }
        }
	init = 1;
    }
    return ct[cc];
}
