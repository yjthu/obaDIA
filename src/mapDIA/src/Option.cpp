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
#include"Option.hpp"


Option::Option(const map<string,string>& opm) :
    d_opm(opm),
    d_rep(true),d_log2(true),d_shared_q(false),
    d_min_f(0),d_max_f(INT_MAX),d_min_q(0),d_level(3),d_nc(0),d_ncolumns(0),d_top_q(INT_MAX),d_dp(INT_MIN),
    d_window(0),d_SDF(NAN),d_min_CV(0),d_min_correl(1),d_minDE(0),d_maxDE(0),d_fudge(-1),d_impute_scale(-1)
{
    istringstream iss(d_opm.find("FILE")->second);
    iss>>d_file;
    if (not ifstream(file().c_str())) throw runtime_error("FILE not found");

    int int0;
    double double0;
    string str0;
    iss.clear(); iss.str(d_opm.find("LEVEL")->second);
    if (iss>>int0) {
        d_level=int0;
        if (level()<1 or level()>3) throw runtime_error("LEVEL?");
    }

    iss.clear(); iss.str(d_opm.find("FUDGE")->second);
    if (iss>>double0) {
        d_fudge=double0;
        if (fudge()<0 or fudge()>1) throw runtime_error("FUDGE?");
    }

    iss.clear(); iss.str(d_opm.find("IMPUTE")->second);
    if (iss>>str0>>double0) {
        transform(str0.begin(),str0.end(),str0.begin(),::tolower);
        d_impute_type=str0;
        if (impute_type()!="row" and impute_type()!="group") throw runtime_error("IMPUTE?");
        d_impute_scale=double0;
        if (impute_scale()<=0) throw runtime_error("IMPUTE?");
    }

    iss.clear(); iss.str(d_opm.find("EXPERIMENTAL_DESIGN")->second);
    iss>>str0;
    transform(str0.begin(),str0.end(),str0.begin(),::tolower);
    if (str0!="independentdesign" and str0!="replicatedesign") {
        throw runtime_error("EXPERIMENTAL_DESIGN?");
    }
    d_rep=(str0=="replicatedesign");

    iss.clear(); iss.str(d_opm.find("LOG2_TRANSFORMATION")->second);
    if (iss>>str0) {
        transform(str0.begin(),str0.end(),str0.begin(),::tolower);
        if (str0!="true" and str0!="false") throw runtime_error("LOG2_TRANSFORMATION");
        d_log2=(str0=="true");
    }

    iss.clear(); iss.str(d_opm.find("REMOVE_SHARED_PEPTIDE")->second);
    if (iss>>str0) {
        transform(str0.begin(),str0.end(),str0.begin(),::tolower);
        if (str0!="true" and str0!="false") throw runtime_error("REMOVE_SHARED_PEPTIDE");
        d_shared_q=(str0=="true");
    }

    iss.clear(); iss.str(d_opm.find("NORMALIZATION")->second);
    iss>>d_normalization;
    transform(normalization().begin(),normalization().end(),d_normalization.begin(),::toupper);
    if (normalization()=="RT") {
        if (not (iss>>d_window) or window()<=0) throw runtime_error("RT window?");
        if (iss>>int0) d_dp=int0;
    }

    if (level()>=2) {
        iss.clear(); iss.str(d_opm.find("SDF")->second);
        if (iss>>double0) {
            d_SDF=double0;
            if (SDF()<=0) throw runtime_error("SDF must be positive");
        }

        iss.clear(); iss.str(d_opm.find("PSEUDOCV")->second);
        if (iss>>double0) {
            d_min_CV=double0;
            if (min_CV()<0 or min_CV()>1) throw runtime_error("PSEUDOCV?");
        }

        iss.clear(); iss.str(d_opm.find("MIN_CORREL")->second);
        if (not (iss>>d_min_correl) or (min_correl()<-1 or min_correl()>1)) throw runtime_error("MIN_CORREL?");
        if (level()==2) {
            iss.clear(); iss.str(d_opm.find("MIN_PEP_PER_PROT")->second);
            if (not (iss>>d_min_f) or min_f()<1) throw runtime_error("MIN_PEP_PER_PROT?");
            iss.clear(); iss.str(d_opm.find("MAX_PEP_PER_PROT")->second);
            if (iss>>int0) {
                d_max_f=int0;
                if (max_f()<min_f()) throw runtime_error("MAX_PEP_PER_PROT < MIN_PEP_PER_PROT");
                d_top_q = d_max_f;
            }
        } else {
            iss.clear(); iss.str(d_opm.find("MIN_FRAG_PER_PEP")->second);
            if (not (iss>>d_min_f) or min_f()<1) throw runtime_error("MIN_FRAG_PER_PEP?");
            iss.clear(); iss.str(d_opm.find("MAX_FRAG_PER_PEP")->second);
            if (iss>>int0) {
                d_max_f=int0;
                if (max_f()<min_f()) throw runtime_error("MAX_FRAG_PER_PEP < MIN_FRAG_PER_PEP");
            }
            iss.clear(); iss.str(d_opm.find("MIN_PEP_PER_PROT")->second);
            if (not (iss>>d_min_q) or min_q()<1) throw runtime_error("MIN_PEP_PER_PROT?");
            iss.clear(); iss.str(d_opm.find("MAX_PEP_PER_PROT")->second);
            if (iss>>int0) {
                d_top_q=int0;
                if (top_q()<min_q()) throw runtime_error("MAX_PEP_PER_PROT < MIN_PEP_PER_PROT");
            }
        }
    }

    iss.clear(); iss.str(d_opm.find("LABELS")->second);
    for (;iss>>str0;d_labels.push_back(str0));

    if (rep()) {
        iss.clear(); iss.str(d_opm.find("MIN_OBS")->second);
        if (not (iss>>int0)) throw runtime_error("MIN_OBS?");
        d_min_obs.assign(labels().size(),int0);
        iss.clear(); iss.str(d_opm.find("SIZE")->second);
        if (not (iss>>int0)) throw runtime_error("SIZE?");
        d_ssize.assign(labels().size(),int0);
    } else {
        iss.clear(); iss.str(d_opm.find("MIN_OBS")->second);
        for (;iss>>int0;d_min_obs.push_back(int0));
        if (labels().size()!=min_obs().size()) throw runtime_error("Unequal dimension for LABELS and MIN_OBS");
        iss.clear(); iss.str(d_opm.find("SIZE")->second);
        for (;iss>>int0;d_ssize.push_back(int0));
        if (labels().size()!=ssize().size()) throw runtime_error("Unequal dimension for LABELS and SIZE");
    }

    for (unsigned i=0;i<min_obs().size();i++) {
        if (min_obs().at(i)>ssize().at(i)) throw runtime_error("MIN_OBS>SIZE");
    }

    iss.clear(); iss.str(d_opm.find("MIN_DE")->second);
    if (not (iss>>d_minDE) or minDE()<0) throw runtime_error("MIN_DE?");
    iss.clear(); iss.str(d_opm.find("MAX_DE")->second);
    if (not (iss>>d_maxDE) or maxDE()<minDE() or maxDE()>1) throw runtime_error("MAX_DE?");


    iss.clear(); iss.str(d_opm.find("INCLUSION_LIST")->second);
    for (;iss>>str0;d_inclusion.insert(str0));

    istringstream iss0(d_opm.find("CONTRAST")->second);
    vector<vector<char> > cindex(nt(),vector<char>(nt()));
    for (int i=0;i<nt();i++) {
        if (not getline(iss0,str0)) throw runtime_error("CONTRAST?");
        iss.clear(); iss.str(str0);
        for (int j=0;j<nt();j++) {
            char c;
            if (iss>>c) cindex.at(i).at(j)=c; else throw runtime_error("CONTRAST?");
        }
    }
    for (int i=0;i<nt();i++) for (int j=0;j<i;j++) {
        if (cindex.at(i).at(j)=='1' and cindex.at(j).at(i)=='1') {
            throw runtime_error("CONTRAST("+to<string>(i+1)+","+to<string>(j+1)+") = CONTRAST("+to<string>(j+1)+","+to<string>(i+1)+") = 1");
        }
    }
    int0=0;
    for (int i=0;i<nt();i++) for (int j=0;j<nt();j++) if (i!=j) {
        if (cindex.at(i).at(j)=='1') {
            int0++;
            d_ct1.push_back(i);
            d_ct2.push_back(j);
        } else if (cindex.at(i).at(j)!='0') {
            throw runtime_error("CONTRAST not binary matrix");
        }
    }
    d_nc=int0;

    d_np.assign(ssize().size(),0);
    partial_sum(ssize().begin(),ssize().end()-1,d_np.begin()+1);
    d_ncolumns=accumulate(d_ssize.begin(),d_ssize.end(),0);
}
