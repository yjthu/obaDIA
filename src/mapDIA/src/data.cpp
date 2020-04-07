/*
   Copyright (C) <2015-2017> Guoshou Teo <guoshou@u.nus.edu> and
   Hyungwon Choi <hyung_won_choi@nuhs.edu.sg>,
   National University of Singapore.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   */


#include"main.hpp"
#include"Module.hpp"
#include"Option.hpp"
#include"Pre.hpp"


bool read_param(ifstream& ifs,string& lstr,string& rstr,const map<string,string>& opm)
{
    istringstream liss;
    for (string str0;getline(ifs,str0);) {
        str0=str0.substr(0,str0.find("#"));
        const size_t e=str0.find("=");
        if (e!=string::npos) {
            liss.str(str0.substr(0,e));
            liss>>lstr;
            if (opm.find(lstr)==opm.end()) throw runtime_error("Unknown parameter: "+lstr);
            if (lstr=="CONTRAST") {
                getline(ifs,rstr);
                istringstream iss(rstr);
                int n=0;
                for (char j;iss>>j;n++);
                for (int i=1;i<n and getline(ifs,str0);i++,rstr+='\n'+str0);
            } else {
                rstr=str0.substr(e+1);
            }
            return true;
        }
    }
    return false;
}


Option set_param(const string& filepath,Module& mo,Module& mo1)
{
    ifstream input_ifs(filepath.c_str());
    if (not input_ifs) throw runtime_error("can't open "+filepath);

    map<string,string> opm;
    opm["FILE"]; opm["LEVEL"]; opm["EXPERIMENTAL_DESIGN"];
    opm["NORMALIZATION"]; opm["SDF"]; opm["MIN_CORREL"];
    opm["MIN_OBS"]; opm["MIN_FRAG_PER_PEP"]; opm["MAX_FRAG_PER_PEP"];
    opm["MIN_PEP_PER_PROT"]; opm["LABELS"]; opm["SIZE"];
    opm["MIN_DE"]; opm["MAX_DE"]; opm["CONTRAST"];
    opm["MAX_PEP_PER_PROT"]; opm["MODULE"]; opm["MODULE_TYPE"];
    opm["MRF_TYPE"]; opm["MODULE_SIZE"]; opm["MODULE_FREQ"];
    opm["MODULE2"]; opm["MODULE_TYPE2"]; opm["MRF_TYPE2"];
    opm["MODULE_SIZE2"]; opm["MODULE_FREQ2"]; opm["LOG2_TRANSFORMATION"];
    opm["FUDGE"]; opm["REMOVE_SHARED_PEPTIDE"]; opm["IMPUTE"];
    opm["INCLUSION_LIST"]; opm["PSEUDOCV"];

    set<string> opset;
    for (string lstr,rstr;read_param(input_ifs,lstr,rstr,opm);opm[lstr]=rstr) {
        if (opset.find(lstr)!=opset.end()) throw runtime_error("Duplicate \""+lstr+" = \"");
        opset.insert(lstr);
    }

    mo.setModulebool(opm.find("MODULE")->second);

    if (mo.modulebool()) {
        mo.setModule(opm.find("MODULE")->second);
        mo.setModule_type(opm.find("MODULE_TYPE")->second);
        mo.setMrf(opm.find("MRF_TYPE")->second);
        if (mo.module_type()=="group_list") {
            mo.setModule_size(opm.find("MODULE_SIZE")->second);
            mo.setModule_freq(opm.find("MODULE_FREQ")->second);
        }

        mo1.setModulebool(opm.find("MODULE2")->second);

        if (mo1.modulebool()) {
            mo1.setModule(opm.find("MODULE2")->second);
            mo1.setModule_type(opm.find("MODULE_TYPE2")->second);
            mo1.setMrf(opm.find("MRF_TYPE2")->second);
            if (mo1.module_type()=="group_list") {
                mo1.setModule_size(opm.find("MODULE_SIZE2")->second);
                mo1.setModule_freq(opm.find("MODULE_FREQ2")->second);
            }
        }
    }

    Option op(opm);
    return op;
}
