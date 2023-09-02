
#include "rtklib.h"

typedef struct {
    double sec;
	double wm[3], vm[3];
} ImuData;
typedef struct {
    double sec;
	int wm[3], vm[3];
} ImuData1;

typedef struct {
    double sec;
	double att, stdatt, vngps[3], stdvn[3], posgps[3], stdpos[3];		
} GnssData;

int main() {
    FILE *fp_in, *fp_out;
    // char infile[1024] = "G:\\mGNSS\\GNSSINS\\emulation\\gnss.bin";
    // char outfile[1024] = "G:\\mGNSS\\GNSSINS\\emulation\\outgnss.txt";
    //char infile[1024] = "G:\\mGNSS\\GNSSINS\\emulation\\ins.bin";
    char infile[1024]  = "/media/zhuhang/D/Data/PPPRTK_INS-Wuh-220519/220520-test/2test/INS-D/F1760369-2022-05-20-13-45-04.imu";
    char outfile[1024] = "/media/zhuhang/D/Data/PPPRTK_INS-Wuh-220519/220520-test/2test/INS-D/INSD-G2.txt";
    char buff[215];
    ImuData imu;
    ImuData1 imu1;
    GnssData gnss;

    fp_in = fopen(infile,"rb");
    fp_out = fopen(outfile,"w");
     while (!feof(fp_in)) {
        if (fread(&imu1, sizeof(ImuData1), 1, fp_in) == 1) {
            // fprintf(fp_out, "%.3f %10d %10d %10d %10d %10d %10d\n", imu1.sec, imu1.wm[0],
            //         imu1.wm[1], imu1.wm[2], imu1.vm[0], imu1.vm[1], imu1.vm[2]);
            fprintf(fp_out, "%10.3f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n", imu1.sec, imu1.wm[0]*1E-04,
                    imu1.wm[1]*1E-04, imu1.wm[2]*1E-04, imu1.vm[0]*1E-04, imu1.vm[1]*1E-04, imu1.vm[2]*1E-04);
            trace(1, "%10.3f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n", imu1.sec, imu1.wm[0]*1E-04,
                    imu1.wm[1]*1E-04, imu1.wm[2]*1E-04, imu1.vm[0]*1E-04, imu1.vm[1]*1E-04, imu1.vm[2]*1E-04);
        }
    }
    fclose(fp_in);
    return 0;

    // fp_in = fopen(infile,"rb");
    // fp_out = fopen(outfile,"w");
    
    // fp_in = fopen("G:\\mGNSS\\GNSSINS\\emulation\\result\\test-imu.txt","rb");
    // fp_out = fopen("G:\\mGNSS\\GNSSINS\\emulation\\result\\result-imu.txt","w");
    // while (!feof(fp_in)) {
    //     if (fread(&imu, sizeof(ImuData), 1, fp_in) == 1) {
    //         fprintf(fp_out, "%.9f %.9f %.9f %.9f %.9f %.9f %.9f\n", imu.sec, imu.wm[0],
    //                 imu.wm[1], imu.wm[2], imu.vm[0], imu.vm[1], imu.vm[2]);
    //     }
    // }

    fp_in = fopen("G:\\mGNSS\\GNSSINS\\emulation\\gnss.bin","rb");
    fp_out = fopen("G:\\mGNSS\\GNSSINS\\emulation\\result\\result-gps.txt","w");
    while (!feof(fp_in)) {
        if (fread(&gnss, sizeof(GnssData), 1, fp_in) == 1) {
            fprintf(fp_out, "%.9f %.9f %.9f %.9f %.9f %.9f %.9f %.9f\n", gnss.att, gnss.vngps[0], 
                    gnss.vngps[1], gnss.vngps[2], gnss.posgps[0], gnss.posgps[1], gnss.posgps[2], gnss.sec);
        }
    }

    fclose(fp_in);
    fclose(fp_out);

    // fp_in = fopen(outfile,"r");
    // fp_out = fopen(result,"wb");

    // while(fgets(buff, sizeof(buff), fp_in)) {
    //     if (sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf",  &imu.sec, imu.wm, 
    //                imu.wm+1, imu.wm+2, imu.vm, imu.vm+1, imu.vm+2) < 7) {
	// 		continue;
	// 	}
    //     fwrite(&imu,sizeof(ImuData), 1, fp_out);
    // }
    // fclose(fp_in);
    // fclose(fp_out);

    trace(1, "success\n");
    getchar();
    return 1;
}