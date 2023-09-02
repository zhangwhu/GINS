
#include "CGNSSConfig.h"

/******************************  class CGNSSConfig  ***********************************/
CGNSSConfig::CGNSSConfig() 
{
    
}

int CGNSSConfig::LoadYAML(const char *file) 
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
            str2time(first_layer->second.as<std::string>().c_str(), 0, 32, &prcopt.ts);
        }
        if (first_layer->first.as<std::string>() == "endtime") {
            str2time(first_layer->second.as<std::string>().c_str(), 0, 32, &prcopt.te);
        }
        if (first_layer->first.as<std::string>() == "soltype") {
            prcopt.soltype = first_layer->second.as<int>();
        }
        if (first_layer->first.as<std::string>() == "simulated-interrupt") {
            
        }
        if (first_layer->first.as<std::string>() == "file") {
            YAML::Node it = first_layer->second.as<YAML::Node>();
            if (it["rovobsfile"])    strcpy(filopt.rovobs,   it["rovobsfile"].as<std::string>().c_str());
            if (it["refobsfile"])    strcpy(filopt.refobs,   it["refobsfile"].as<std::string>().c_str());

            if (it["satantfile"])    strcpy(filopt.satantp,  it["satantfile"].as<std::string>().c_str());
            if (it["rcvantfile"])    strcpy(filopt.rcvantp,  it["rcvantfile"].as<std::string>().c_str());
            if (it["staposfile"])    strcpy(filopt.stapos,   it["staposfile"].as<std::string>().c_str());
            if (it["brdcfile"])      strcpy(filopt.brdc,     it["brdcfile"].as<std::string>().c_str());
            if (it["clkfile"])       strcpy(filopt.clk,      it["clkfile"].as<std::string>().c_str());
            if (it["sp3file"])       strcpy(filopt.sp3,      it["sp3file"].as<std::string>().c_str());
            if (it["dcbfile"])       strcpy(filopt.dcb,      it["dcbfile"].as<std::string>().c_str());
            if (it["eopfile"])       strcpy(filopt.eop,      it["eopfile"].as<std::string>().c_str());
            if (it["blqfile"])       strcpy(filopt.blq,      it["blqfile"].as<std::string>().c_str());
            if (it["fcbfile"])       strcpy(filopt.fcb,      it["fcbfile"].as<std::string>().c_str());
            if (it["corrfile"])      strcpy(filopt.corr,     it["corrfile"].as<std::string>().c_str());

            if (it["outpath"])       strcpy(filopt.outdir,   it["outpath"].as<std::string>().c_str());
            if (it["outfile1"])      strcpy(filopt.outfile1, it["outfile1"].as<std::string>().c_str());
            if (it["outfile2"])      strcpy(filopt.outfile2, it["outfile2"].as<std::string>().c_str());
            if (it["tracefile"])     strcpy(filopt.trace,    it["tracefile"].as<std::string>().c_str());
            if (it["geexefile"])     strcpy(filopt.geexe,    it["geexefile"].as<std::string>().c_str());
        }

        if (first_layer->first.as<std::string>() == "gnss") {
            char exsats_[1024];
            YAML::Node it = first_layer->second.as<YAML::Node>();

            if (it["pos1-posmode"])      prcopt.mode     =   it["pos1-posmode"].as<int>();
            if (it["pos1-interval"])     prcopt.ti       =   it["pos1-interval"].as<double>();
            if (it["pos1-reinit"])       prcopt.reinit   =   it["pos1-reinit"].as<double>();
            if (it["pos1-frequency"])    prcopt.nf       =   it["pos1-frequency"].as<int>();
            if (it["pos1-elmask"])       prcopt.elmin    =   it["pos1-elmask"].as<double>()*D2R;
            if (it["pos1-snrmask_r"])    prcopt.snrmask.ena[0]    =     it["pos1-snrmask_r"].as<bool>();
            if (it["pos1-snrmask_b"])    prcopt.snrmask.ena[1]    =     it["pos1-snrmask_b"].as<bool>();
            if (it["pos1-snrmask_L1"])   matcpy(prcopt.snrmask.mask[0], it["pos1-snrmask_L1"].as<std::vector<double>>().data(),9,1);
            if (it["pos1-snrmask_L2"])   matcpy(prcopt.snrmask.mask[1], it["pos1-snrmask_L2"].as<std::vector<double>>().data(),9,1);
            if (it["pos1-snrmask_L3"])   matcpy(prcopt.snrmask.mask[2], it["pos1-snrmask_L2"].as<std::vector<double>>().data(),9,1);
            if (it["pos1-dynamics"])     prcopt.dynamics =   it["pos1-dynamics"].as<int>();
            if (it["pos1-tidecorr"])     prcopt.tidecorr =   it["pos1-tidecorr"].as<int>();
            if (it["pos1-ionoopt"])      prcopt.ionoopt  =   it["pos1-ionoopt"].as<int>();
            if (it["pos1-tropopt"])      prcopt.tropopt  =   it["pos1-tropopt"].as<int>();
            if (it["pos1-sateph"])       prcopt.sateph   =   it["pos1-sateph"].as<int>();
            if (it["pos1-posopt"])       intcpy(prcopt.posopt,it["pos1-posopt"].as<std::vector<int>>().data(),9,1);
            if (it["pos1-exclsats"])     strcpy(exsats_,     it["pos1-exclsats"].as<std::string>().c_str());
            if (it["pos1-navsys"])       prcopt.navsys   =   it["pos1-navsys"].as<int>();
            if (it["pos1-outopt"])       prcopt.outopt   =   it["pos1-outopt"].as<int>();

            if (it["pos2-maxage"])       prcopt.maxtdiff =   it["pos2-maxage"].as<double>();
            if (it["pos2-rejion"])       prcopt.maxrej[1]=   it["pos2-rejion"].as<int>();
            if (it["pos2-rejifb"])       prcopt.maxrej[2]=   it["pos2-rejifb"].as<int>();
            if (it["pos2-rejamb"])       prcopt.maxrej[3]=   it["pos2-rejamb"].as<int>();
            if (it["pos2-acctime"])      prcopt.acctime  =   it["pos2-acctime"].as<int>();
            if (it["pos2-oution"])       prcopt.maxout[1]=   it["pos2-oution"].as<int>();
            if (it["pos2-outifb"])       prcopt.maxout[2]=   it["pos2-outifb"].as<int>();
            if (it["pos2-outamb"])       prcopt.maxout[3]=   it["pos2-outamb"].as<int>();
            if (it["pos2-cslipgf12"])    prcopt.thresslipgf[0] = it["pos2-cslipgf12"].as<double>();
            if (it["pos2-cslipmw12"])    prcopt.thresslipmw[0] = it["pos2-cslipmw12"].as<double>();
            if (it["pos2-cslipgf23"])    prcopt.thresslipgf[1] = it["pos2-cslipgf23"].as<double>();
            if (it["pos2-cslipmw23"])    prcopt.thresslipmw[1] = it["pos2-cslipmw23"].as<double>();
            if (it["pos2-maxgdop"])      prcopt.maxgdop  =   it["pos2-maxgdop"].as<double>();
            if (it["pos2-prithres"])     prcopt.threscheck[0]  = it["pos2-prithres"].as<double>();
            if (it["pos2-prothres"])     prcopt.threscheck[1]  = it["pos2-prothres"].as<double>();
            if (it["pos2-tdcthres"])     prcopt.threscheck[2]  = it["pos2-tdcthres"].as<double>();
            if (it["pos2-p12thres"])     prcopt.threscheck[3]  = it["pos2-p12thres"].as<double>();
            if (it["pos2-sppthres"])     prcopt.threscheck[4]  = it["pos2-sppthres"].as<double>();
            if (it["pos2-arthres"])      prcopt.threscheck[5]  = it["pos2-arthres"].as<double>();
            if (it["pos2-baselinelen"])  prcopt.baseline[0]    = it["pos2-baselinelen"].as<double>();
            if (it["pos2-baselinesig"])  prcopt.baseline[1]    = it["pos2-baselinesig"].as<double>();

            if (it["pos3-armode"])       prcopt.modear   =   it["pos3-armode"].as<int>();
            if (it["pos3-artype"])       prcopt.typear   =   it["pos3-artype"].as<int>();
            if (it["pos3-arelmask"])     prcopt.elmaskar =   it["pos3-arelmask"].as<double>()*D2R;
            if (it["pos3-arlockcnt"])    prcopt.minlock  =   it["pos3-arlockcnt"].as<int>();
            if (it["pos3-arminsats"])    prcopt.minfixsats=  it["pos3-arminsats"].as<int>();
            if (it["pos3-arthresel"])    prcopt.thresar[0]=  it["pos3-arthresel"].as<double>();
            if (it["pos3-arthreswl"])    prcopt.thresar[1]=  it["pos3-arthreswl"].as<double>();
            if (it["pos3-arthresnl"])    prcopt.thresar[2]=  it["pos3-arthresnl"].as<double>();
            if (it["pos3-arratio"])      prcopt.thresar[3]=  it["pos3-arratio"].as<double>();

            if (it["ant1-postype"])      prcopt.rovpos   =   it["ant1-postype"].as<int>();
            if (it["ant1-pos1"])         prcopt.ru[0]    =   it["ant1-pos1"].as<double>();
            if (it["ant1-pos2"])         prcopt.ru[1]    =   it["ant1-pos2"].as<double>();
            if (it["ant1-pos3"])         prcopt.ru[2]    =   it["ant1-pos3"].as<double>();
            if (it["ant1-anttype"])      strcpy(prcopt.anttype[0],it["ant1-anttype"].as<std::string>().c_str());
            if (it["ant1-antdele"])      prcopt.antdel[0][0]  =   it["ant1-antdele"].as<double>();
            if (it["ant1-antdeln"])      prcopt.antdel[0][1]  =   it["ant1-antdeln"].as<double>();
            if (it["ant1-antdelu"])      prcopt.antdel[0][2]  =   it["ant1-antdelu"].as<double>();

            if (it["ant2-postype"])      prcopt.refpos   =   it["ant2-postype"].as<int>();
            if (it["ant2-pos1"])         prcopt.rb[0]    =   it["ant2-pos1"].as<double>();
            if (it["ant2-pos2"])         prcopt.rb[1]    =   it["ant2-pos2"].as<double>();
            if (it["ant2-pos3"])         prcopt.rb[2]    =   it["ant2-pos3"].as<double>();
            if (it["ant2-anttype"])      strcpy(prcopt.anttype[1],it["ant2-anttype"].as<std::string>().c_str());
            if (it["ant2-antdele"])      prcopt.antdel[1][0]  =   it["ant2-antdele"].as<double>();
            if (it["ant2-antdeln"])      prcopt.antdel[1][1]  =   it["ant2-antdeln"].as<double>();
            if (it["ant2-antdelu"])      prcopt.antdel[1][2]  =   it["ant2-antdelu"].as<double>();

            if (it["misc-timeinterp"])   prcopt.intpref  =        it["misc-timeinterp"].as<bool>();
            if (it["misc-rnxopt1"])      strcpy(prcopt.rnxopt[0], it["misc-rnxopt1"].as<std::string>().c_str());
            if (it["misc-rnxopt2"])      strcpy(prcopt.rnxopt[1], it["misc-rnxopt2"].as<std::string>().c_str());

            if (it["out-solformat"])     solopt.posf     =   it["out-solformat"].as<int>();
            if (it["out-outhead"])       solopt.outhead  =   it["out-outhead"].as<bool>();
            if (it["out-outopt"])        solopt.outopt   =   it["out-outopt"].as<bool>();
            if (it["out-outvel"])        solopt.outvel   =   it["out-outvel"].as<bool>();
            if (it["out-timesys"])       solopt.times    =   it["out-timesys"].as<int>();
            if (it["out-timeform"])      solopt.timef    =   it["out-timeform"].as<int>();
            if (it["out-timendec"])      solopt.timeu    =   it["out-timendec"].as<int>();
            if (it["out-degform"])       solopt.degf     =   it["out-degform"].as<int>();
            if (it["out-fieldsep"])      strcpy(solopt.sep,  it["out-fieldsep"].as<std::string>().c_str());
            if (it["out-maxsolstd"])     solopt.maxsolstd =  it["out-maxsolstd"].as<double>();
            if (it["out-height"])        solopt.height   =   it["out-height"].as<int>();
            if (it["out-geoid"])         solopt.geoid    =   it["out-geoid"].as<int>();
            if (it["out-solstatic"])     solopt.solstatic=   it["out-solstatic"].as<int>();
            if (it["out-nmeaintv1"])     solopt.nmeaintv[0]= it["out-nmeaintv1"].as<double>();
            if (it["out-nmeaintv2"])     solopt.nmeaintv[1]= it["out-nmeaintv2"].as<double>();
            if (it["out-outstat"])       solopt.sstat    =   it["out-outstat"].as<int>();

            /* excluded satellites */
            char buff[1024], *p, *id;
            int i, sat;
            for (i=0;i<MAXSAT;i++) prcopt.exsats[i]=0;
                if (exsats_[0]!='\0') {
                strcpy(buff, exsats_);
                for (p=strtok(buff," "); p; p=strtok(NULL," ")) {
                    if (*p=='+') id=p+1; else id=p;
                    if (!(sat=satid2no(id))) continue;
                    prcopt.exsats[sat-1]=*p=='+'?2:1;
                }
            }
        }

        if (first_layer->first.as<std::string>() == "gnss-stochastic") {
            YAML::Node it = first_layer->second.as<YAML::Node>();
            if (it["P0-pos"])            prcopt.P0[0] = it["P0-pos"].as<double>();
            if (it["P0-clk"])            prcopt.P0[1] = it["P0-clk"].as<double>();
            if (it["P0-trp"])            prcopt.P0[2] = it["P0-trp"].as<double>();
            if (it["P0-ion"])            prcopt.P0[3] = it["P0-ion"].as<double>();
            if (it["P0-ifb"])            prcopt.P0[4] = it["P0-ifb"].as<double>();
            if (it["P0-amb"])            prcopt.P0[5] = it["P0-amb"].as<double>();
            if (it["P0-vel"])            prcopt.P0[6] = it["P0-vel"].as<double>();
            if (it["P0-acc"])            prcopt.P0[7] = it["P0-acc"].as<double>();

            if (it["Qt-pos"])            prcopt.Qt[0] = it["Qt-pos"].as<double>();
            if (it["Qt-trp"])            prcopt.Qt[2] = it["Qt-trp"].as<double>();
            if (it["Qt-ion"])            prcopt.Qt[3] = it["Qt-ion"].as<double>();
            if (it["Qt-ifb"])            prcopt.Qt[4] = it["Qt-ifb"].as<double>();
            if (it["Qt-amb"])            prcopt.Qt[5] = it["Qt-amb"].as<double>();
            if (it["Qt-acch"])           prcopt.Qt[6] = it["Qt-acch"].as<double>();
            if (it["Qt-accv"])           prcopt.Qt[7] = it["Qt-accv"].as<double>();

            if (it["Rt-eratio"])         matcpy(prcopt.eratio, it["Rt-eratio"].as<std::vector<double>>().data(), 3, 1);
            if (it["Rt-phase"])          prcopt.err[1] = it["Rt-phase"].as<double>();
            if (it["Rt-phaseel"])        prcopt.err[2] = it["Rt-phaseel"].as<double>();
            if (it["Rt-phasebl"])        prcopt.err[3] = it["Rt-phasebl"].as<double>();
            if (it["Rt-doppler"])        prcopt.err[4] = it["Rt-doppler"].as<double>();
            if (it["Rt-trop"])           prcopt.err[5] = it["Rt-trop"].as<double>();
            if (it["Rt-iono"])           prcopt.err[6] = it["Rt-iono"].as<double>();
        }
    }
    strcpy(sta.name,"Sta1");
    strcpy(sta.antdes, prcopt.anttype[0]);
    sta.deltype = 0;
    matcpy(sta.del, prcopt.antdel[0], 3, 1);
    return 1; 
}