/*
 *      Module:         int_convert
 *      Author:         Freddie Akeroyd, ISIS
 *      Purpose:        Routines to convert from VAX to local integer representations
 *
 */
#include "geniedefs.h"
#include "int_convert.h"

/*
 * Byte swaps for int and short
 */

static const char* rcsid = "@(#)$Id: int_convert.c,v 1.5 1999/01/07 13:10:39 faa Exp $";
USE_RCSID(rcsid);

#if 0
inline unsigned swap_int(unsigned a)
{
    union { unsigned u; unsigned char c[4]; } temp;
    unsigned char ctemp;
    temp.u = a;
    ctemp = temp.c[0]; temp.c[0] = temp.c[3]; temp.c[3] = ctemp; 
    ctemp = temp.c[1]; temp.c[1] = temp.c[2]; temp.c[2] = ctemp; 
    return temp.u;
}
#endif

#define swap_int(a) ( ((a) << 24) | \
                      (((a) << 8) & 0x00ff0000) | \
                      (((a) >> 8) & 0x0000ff00) | \
        	      ((unsigned long)(a) >>24) )

#define swap_short(a) ( ((a & 0xff) << 8) | ((unsigned short)(a) >> 8) )

/* VAXes are little endian */

unsigned short GENIEMECHANISM local_to_vax_short(const unsigned short* s)
{
#if defined(WORDS_BIGENDIAN)
    return swap_short(*s);
#else 
    return *s;
#endif /* WORDS_BIGENDIAN */
}

unsigned short GENIEMECHANISM vax_to_local_short(const unsigned short* s)
{
#if defined(WORDS_BIGENDIAN)
    return swap_short(*s);
#else
    return *s;
#endif /* WORDS_BIGENDIAN */
}

unsigned GENIEMECHANISM local_to_vax_int(const unsigned* i)
{
#if defined(WORDS_BIGENDIAN)
    return swap_int(*i);
#else
    return *i;
#endif /* WORDS_BIGENDIAN */
}

unsigned GENIEMECHANISM vax_to_local_int(const unsigned* i)
{
#if defined(WORDS_BIGENDIAN)
    return swap_int(*i);
#else
    return *i;
#endif /* WORDS_BIGENDIAN */
}

void GENIEMECHANISM local_to_vax_shorts(unsigned short* sa, const int* n)
{
#if defined(WORDS_BIGENDIAN)
    int i;
    for(i=0; i<*n; i++)
    {
	sa[i] = swap_short(sa[i]);
    }
#endif /* WORDS_BIGENDIAN */
    return;
}

void GENIEMECHANISM vax_to_local_shorts(unsigned short* sa, const int* n)
{
#if defined(WORDS_BIGENDIAN)
    int i;
    for(i=0; i<*n; i++)
    {
	sa[i] = swap_short(sa[i]);
    }
#endif /* WORDS_BIGENDIAN */
    return;
}

void GENIEMECHANISM local_to_vax_ints(unsigned* ia, const int* n)
{
#if defined(WORDS_BIGENDIAN)
    int i;
    for(i=0; i<*n; i++)
    {
	ia[i] = swap_int(ia[i]);
    }
#endif /* WORDS_BIGENDIAN */
    return;
}

void GENIEMECHANISM vax_to_local_ints(unsigned* ia, const int* n)
{
#if defined(WORDS_BIGENDIAN)
    int i;
    for(i=0; i<*n; i++)
    {
	ia[i] = swap_int(ia[i]);
    }
#endif /* WORDS_BIGENDIAN */
    return;
}

