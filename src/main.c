
#include <stdarg.h>
#include "rtklib.h"
//const char stime[32] = "2021 08 11 04 00 00.0";
//const char etime[32] = "2021 08 12 00 00 00.0";
const char stime[32] = "2022 02 02 00 00 00.0";
const char etime[32] = "2022 02 02 00 00 00.0";
/* post-process main ------------------------------------------------------------ */
int main(int argc, char **argv)									 
{
	int ret = 0;
	char conf[1024] = "/media/zhuhang/ZH/mGNSS/ppp_static-linux.conf";
	//char conf[1024] = "G:\\mGNSS\\ppp_static.conf";
	gtime_t ts, te;
	prcopt_t prcopt = prcopt_default;																				  
	solopt_t solopt = solopt_default;		  
	filopt_t filopt = { "" };	

	if (!loadopts(conf, sysopts)) return -1;			
	getsysopts(&prcopt, &solopt, &filopt);

	str2time(stime, 0, 32, &ts);
	str2time(etime, 0, 32, &te);

	//ret = loadobs(ts, te, &filopt);

	ret = postpos(&prcopt, &solopt, &filopt);

	//ret = fcbpos(ts, te, &filopt);		

	//ret = ionoModel(&filopt);			

	if (!ret) fprintf(stderr, "%40s\r", "");
	printf("success\n");
	getchar();
	return ret;
}