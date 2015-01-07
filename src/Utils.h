#ifndef UTILS_H
#define UTILS_H

#include <assert.h>

/**
 * Inspired by : http://stackoverflow.com/a/3767904
 *
 */

#ifndef NDEBUG
#define ASSERT_MESSAGE(condition, message)\
	if( !( condition ) ) { fprintf(stderr, message); }\
	assert(condition)
#else
#define ASSERT_MESSAGE(condition, message)\
		((void)0)
#endif
	
#endif
