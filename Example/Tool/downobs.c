
#include "rtklib.h"

int main(int argc, char *argv[]) 
{
    prcopt_t prcopt = prcopt_default;																  
	solopt_t solopt = solopt_default;		  
	filopt_t filopt = { "" };

    /* load conf */
	if (!loadopts(argv[1], sysopts)) return -1;			
	getsysopts(&prcopt, &solopt, &filopt);

    loadobs(prcopt.ts, prcopt.te, filopt.stapos, filopt.outdir);

    getchar();
}