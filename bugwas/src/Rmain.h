#ifndef _RMAIN_H_
#define _RMAIN_H_
#include <R.h>

extern "C" {
	void Rmain(int* argv_nrow, int* argv_ncol, char** argv, int* EXIT_STATUS);
}

#endif

