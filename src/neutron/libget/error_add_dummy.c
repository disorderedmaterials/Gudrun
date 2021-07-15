/*
 * Place holders for Open GENIE error handling routines
 *
 * Freddie Akeroyd, CCLRC ISIS facility, 7/5/98
 */
#include "geniedefs.h"
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

static const char* rcsid = "@(#)$Id: error_add_dummy.c,v 1.3 1999/02/04 13:37:39 faa Exp $";
USE_RCSID(rcsid);

void error_add(const char *obj, const char *fail, const char *soln)
{
	if (obj == 0 || *obj == '\0')
	    return;

	if (strcasecmp(obj, "INFORMATION") != 0)
	{
	    printf("ERROR returned from: %s\n%s\n", obj, fail);
	}
	else
	{
	    printf("%s\n", fail);
	}
        if (soln != 0 && *soln != '\0')
	    printf("SOLUTION: %s\n", soln);
}

void g_printef(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

void gprintdf(const char *fmt, ...)
{
    return;
}

void gprintif(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}
