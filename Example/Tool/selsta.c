
#include "rtklib.h"

#define MAXRNXLEN   (16*MAXOBSTYPE+4) 

static obs_t obss;

int main() 
{
    gtime_t ts = { 0 }, te = { 0 };
    sta_t sta;
    int i, k, n;
    char infile[1024], outfile[1024], *files[MAXEXFILE] = { 0 }, name[8];

    sprintf( infile, "%s%s", "D:\\Data\\PPPRTK_INS-Wuh-220519\\MGEX-FCB\\","*140*");
    sprintf(outfile, "%s%s", "D:\\Data\\PPPRTK_INS-Wuh-220519\\","station.stas");
    for (i = 0; i < MAXEXFILE; i++) {
		if (!(files[i] = (char *)malloc(1024))) {
			for (i--; i >= 0; i--) free(files[i]);
			return 0;
		}
	}
	n = expath(infile, files, MAXEXFILE);

    FILE *fp = fopen(outfile,"w");
    for (i = 0; i < n; i++) {
        int stat = readrnx(files[i], 0, " ", &obss, NULL, &sta);
        for (k=0; k<1024 && k<obss.n; k++) {
            if ((satsys(obss.data[k].sat, NULL)&SYS_CMP)) {
                if(obss.data[k].P[2]!=0.0) {
                    strncpy(name, sta.name, 4);
                    fprintf(fp, "%s %15.4f %15.4f %15.4f\n", name, sta.pos[0], sta.pos[1], sta.pos[2]);
                    break;
                }
            }
        }
        fflush(fp);
        freeobs(&obss);
    }
    fclose(fp);

    for (i = 0; i < MAXEXFILE; i++) free(files[i]);
}