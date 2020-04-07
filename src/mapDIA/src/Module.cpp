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


void Module::setModulebool(const string& str0)
{
    istringstream iss(str0);
    string str1;
    iss>>str1;
    d_modulebool=ifstream(str1.c_str());
}


void Module::setModule_freq(const string& str0)
{
    istringstream iss(str0);
    if (not (iss>>d_module_freq)) throw runtime_error("MODULE_FREQ?");
}


void Module::setModule_size(const string& str0)
{
    istringstream iss(str0);
    if (not (iss>>d_min_size)) throw runtime_error("MODULE_SIZE?");
    int int0;
    if (iss>>int0) d_max_size=int0;
}


void Module::setModule(const string& str0)
{
    istringstream iss(str0);
    iss>>d_module;
}


void Module::setModule_type(const string& str0)
{
    istringstream iss(str0);
    iss>>d_module_type;
    transform(module_type().begin(),module_type().end(),d_module_type.begin(),::tolower);
    if (module_type()!="group_list" and module_type()!="edge_list") throw runtime_error("MODULE_TYPE?");
}


void Module::setMrf(const string& str0)
{
    istringstream iss(str0);
    iss>>d_mrf;
    if (mrf()!="0_1" and mrf()!="-1_1") throw runtime_error("MRF_TYPE?");
}


void Module::module_data(const vector<string>& pidvec)
{
    ifstream module_ifs(module().c_str());
    if (not module_ifs) throw runtime_error("can't open "+module());
    if (module_type()=="edge_list") {
        d_adj.resize(pidvec.size());
        for (string str0,str1;getline(module_ifs,str0);) {
            istringstream iss(str0);
            if (getline(iss,str0,'\t') and getline(iss,str1,'\t') and str0!=str1) {
                const unsigned i=find(pidvec.begin(),pidvec.end(),str0)-pidvec.begin();
                const unsigned j=find(pidvec.begin(),pidvec.end(),str1)-pidvec.begin();
                if (i!=pidvec.size() and j!=pidvec.size()) {//if both found
                    d_adj.at(j).insert(i);
                    d_adj.at(i).insert(j);
                }
            }
        }
    }
    if (module_type()=="group_list") {
        d_adj.resize(pidvec.size());
        multimap<string,string> GO_mm;
        for (string str0,str1;getline(module_ifs,str0);) if (str0.at(0)!='!') {
            istringstream iss(str0);
            getline(iss,str0,'\t');//skip DB
            getline(iss,str1,'\t');//read DB Object ID
            getline(iss,str0,'\t');//skip DB Object Symbol
            getline(iss,str0,'\t');//skip Qualifier
            getline(iss,str0,'\t');//read GO ID
            GO_mm.insert(make_pair(str0,str1));
        }
        vector<vector<int> > GOfreq(pidvec.size());
        for (unsigned i=0;i<pidvec.size();i++) GOfreq.at(i).assign(i,0);
        for (multimap<string,string>::const_iterator end1,itg=GO_mm.begin();itg!=GO_mm.end();itg=end1) {
            end1=GO_mm.upper_bound(itg->first);
            int GOsize=0;
            for (multimap<string,string>::const_iterator it=itg;it!=end1;it++,GOsize++);
            if (GOsize<min_size() or GOsize>max_size()) continue;
            set<int> GOset;
            for (multimap<string,string>::const_iterator it=itg;it!=end1;it++) {
                const unsigned i=find(pidvec.begin(),pidvec.end(),it->second)-pidvec.begin();
                if (i!=pidvec.size()) GOset.insert(i);
            }
            if (GOset.size()>1) {
                set<int>::const_iterator it0=GOset.begin(),it1;
                for (it0++;it0!=GOset.end();it0++) for (it1=GOset.begin();it1!=it0;it1++) GOfreq.at(*it0).at(*it1)++;
            }
        }

        for (unsigned i=1;i<GOfreq.size();i++) for (unsigned j=0;j<i;j++) if (GOfreq.at(i).at(j)>=module_freq()) {
            d_adj.at(j).insert(i);
            d_adj.at(i).insert(j);
        }
    }

    int sum=0;
    for (unsigned p=0;p<pidvec.size();p++) sum += adj().at(p).size();
    cout<<sum/2.<<" edges\n";
}
