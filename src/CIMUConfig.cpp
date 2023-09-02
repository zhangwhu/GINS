
#include "CIMUConfig.h"

/******************************  class CIMUConfig  ***********************************/
CIMUConfig::CIMUConfig() 
{
    imuopt = { 0 };
}

static double epoch2gpst(const double *ep, int *week)
{
    double time, frac, dt;
	const int doy[] = { 1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 };
    const int time0 = 315964800;
	int days, sec, year = (int)ep[0], mon = (int)ep[1], day = (int)ep[2];

	if (year<1970 || 2099<year || mon<1 || 12<mon) return 0;

	/* leap year if year%4==0 in 1901-2099 */
	days = (year - 1970) * 365 + (year - 1969) / 4 + doy[mon - 1] + day - 2 + (year % 4 == 0 && mon >= 3 ? 1 : 0);
	sec  = (int)floor(ep[5]);
	time = days * 86400 + (int)ep[3] * 3600 + (int)ep[4] * 60 + sec;
	frac  = ep[5] - sec;

    dt = time -time0;
    int w = (int)(dt / (86400 * 7));
    if (week) *week = w;
	return (double)(dt - (double)w * 86400 * 7) + frac;
}

static int str2gpst(const char *s, int i, int n, int *week, double *tow)
{
	double ep[6];
	char str[256], *p = str;
	int mon, day, sec;
	const int mday[] = { /* # of days in a month */
		31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31,
		31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
	};

	if (i<0 || (int)strlen(s)<i || (int)sizeof(str) - 1<i) return -1;
	for (s += i; *s&&--n >= 0;) *p++ = *s++;
	*p = '\0';
	if (sscanf(str, "%lf %lf %lf %lf %lf %lf", ep, ep + 1, ep + 2, ep + 3, ep + 4, ep + 5) == 6 ||
		sscanf(str, "%lf/%lf/%lf %lf:%lf:%lf", ep, ep + 1, ep + 2, ep + 3, ep + 4, ep + 5) == 6){
		if (ep[0]<100.0) ep[0] += ep[0]<80.0 ? 2000.0 : 1900.0;
	}
	else if (sscanf(str, "%lf:%lf:%lf", ep, ep + 1, ep + 3) == 3) {
		if (ep[0]<100.0) ep[0] += ep[0]<80.0 ? 2000.0 : 1900.0;
		day = (int)ep[1];
		for (mon = 0; mon<48; mon++) {
			if (day >= mday[mon]) day -= mday[mon]; else break;
		}
		ep[1] = mon % 12 + 1; ep[2] = day;
		sec = (int)ep[3];
		ep[3] = sec / 3600; ep[4] = sec % 3600 / 60; ep[5] = sec % 60;
	}
	else return -1;

	*tow = epoch2gpst(ep, week);
	return 0;
}

int CIMUConfig::LoadYAML(const char *file) 
{
	YAML::Node config;
	try {
		config = YAML::LoadFile(file);
	}
	catch (YAML::Exception &exception) {
		printf("Failed to read configuration file\n");
		return 0;
	}
    
    for (YAML::const_iterator first_layer = config.begin(); first_layer != config.end(); ++first_layer) {
        if (first_layer->first.as<std::string>() == "statime") {
            str2gpst(first_layer->second.as<std::string>().c_str(), 0, 32, NULL, &imuopt.ts);
        }
        if (first_layer->first.as<std::string>() == "endtime") {
            str2gpst(first_layer->second.as<std::string>().c_str(), 0, 32, NULL, &imuopt.te);
        }
        if (first_layer->first.as<std::string>() == "imu2gps") {
            imuopt.imu2gps = first_layer->second.as<int>();
        }
        if (first_layer->first.as<std::string>() == "simulated-interrupt") {
            
        }
        if (first_layer->first.as<std::string>() == "file") {
            YAML::Node it = first_layer->second.as<YAML::Node>();
            if (it["imufile"])       strcpy(imuopt.imufile,  it["imufile"].as<std::string>().c_str());
            if (it["outpath"])       strcpy(imuopt.outpath,  it["outpath"].as<std::string>().c_str());
        }

        if (first_layer->first.as<std::string>() == "imu") {
            YAML::Node it = first_layer->second.as<YAML::Node>();
            if (it["data-freq"])         imuopt.freq      =  it["data-freq"].as<int>();
            if (it["data-rfu"])          strncpy(imuopt.rfu, it["data-rfu"].as<std::string>().c_str(),3);

            if (it["gyr-unit"])          imuopt.gyr_unit  =  it["gyr-unit"].as<int>();
            if (it["gyr-type"])          imuopt.gyr_type  =  it["gyr-type"].as<int>();
            if (it["acc-unit"])          imuopt.acc_unit  =  it["acc-unit"].as<int>();
            if (it["acc-type"])          imuopt.acc_type  =  it["acc-type"].as<int>();

            if (it["aln-mode"])          imuopt.aln_mode  =  it["aln-mode"].as<int>();
            if (it["aln-t1"])            imuopt.aln_t1    =  it["aln-t1"].as<double>();
            if (it["aln-t2"])            imuopt.aln_t2    =  it["aln-t2"].as<double>();
        }

        if (first_layer->first.as<std::string>() == "installation") {
            YAML::Node it = first_layer->second.as<YAML::Node>();
            if (it["antlever"])          imuopt.antlever =   it["antlever"].as<std::vector<double>>().data();
            if (it["bodyangle"])         imuopt.bodyangle=   it["bodyangle"].as<std::vector<double>>().data();
        }

        if (first_layer->first.as<std::string>() == "imu-stochastic") {
            YAML::Node it = first_layer->second.as<YAML::Node>();
            CMat3 pUnit(1.0/glv.Re, 1.0/glv.Re, 1.0);
            if (it["P0-phi"])            imuopt.P0_phi = CVect3(it["P0-phi"].as<std::vector<double>>().data())*glv.deg;
            if (it["P0-gyr"])            imuopt.P0_gyr = CVect3(it["P0-gyr"].as<std::vector<double>>().data())*glv.dps;
            if (it["P0-acc"])            imuopt.P0_acc = it["P0-acc"].as<std::vector<double>>().data();
            if (it["P0-vel"])            imuopt.P0_vel = it["P0-vel"].as<std::vector<double>>().data();
            if (it["P0-pos"])            imuopt.P0_pos = CVect3(it["P0-pos"].as<std::vector<double>>().data())*pUnit;

            if (it["Qt-phi"])            imuopt.Qt_phi = CVect3(it["Qt-phi"].as<std::vector<double>>().data())*glv.deg;
            if (it["Qt-gyr"])            imuopt.Qt_gyr = CVect3(it["Qt-gyr"].as<std::vector<double>>().data())*glv.dps;
            if (it["Qt-acc"])            imuopt.Qt_acc = it["Qt-acc"].as<std::vector<double>>().data();
            if (it["Qt-vel"])            imuopt.Qt_vel = it["Qt-vel"].as<std::vector<double>>().data();
            if (it["Qt-pos"])            imuopt.Qt_pos = CVect3(it["Qt-pos"].as<std::vector<double>>().data())*pUnit;

            // if (it["P0-lv"])             imuopt.P0_lv  = it["P0-lv"].as<double>();
            // if (it["P0-KG"])             imuopt.P0_KG  = it["P0-KG"].as<std::vector<double>>().data();
            // if (it["P0-KA"])             imuopt.P0_KA  = it["P0-KA"].as<std::vector<double>>().data();
        }

        if (first_layer->first.as<std::string>() == "filter") {
            YAML::Node it = first_layer->second.as<YAML::Node>();
            if (it["nSample"])           imuopt.nn     = it["nSample"].as<int>();
            if (it["TDKF-Steps"])        imuopt.tdStep = it["TDKF-Steps"].as<int>();
            if (it["nStates"])           imuopt.nq     = it["nStates"].as<int>();
        }
    }
    return 1; 
}