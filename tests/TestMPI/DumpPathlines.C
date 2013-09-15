#include <stdio.h>
#include "OSUFlow.h"

#include "PathlineLoader.h"

int main(int argc, char **argv)
{
	PathlineLoader trace(argv[1]);
	trace.connectTraces();
	trace.dump();
	return 0;
}
