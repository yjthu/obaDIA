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
#include"Post.hpp"
#include"Est.hpp"


Option set_param(const string& filepath,Module& mo,Module& mo1);

void print_data(const Pre& pr);
void print_selection(const Pre& pr);
void print_analysis(const Pre& pr,const Est& es,const Est& es_);


int main(int argc, char** argv)
{
    if (argc==1 or string(argv[ 1 ])=="--help") {
        cout<<"Usage: "+string(argv[ 0 ])+" [INPUT PARAMETER FILE]\n";
        return 0;
    }
    if (string(argv[ 1 ])=="--version") {
        cout<<"mapDIA 3.1.0\n";
        return 0;
    }

    Module mo,mo1;
    const Option& op=set_param(string(argv[ 1 ]),mo,mo1);

    Pre pr(op);

    print_selection(pr);


    print_data(pr);

    if (op.nc()>0) {
        Post po(pr);

        Est es(po);

        if (mo.modulebool()) {
            mo.module_data(pr.pidvec());
            if (mo1.modulebool()) mo1.module_data(pr.pidvec());

            vector<int> pbool(po.yy().size());
            for (unsigned p=0;p<po.yy().size();p++) {
                if (mo.adj().at(p).size()>0 or (mo1.modulebool() and mo1.adj().at(p).size()>0)) {
                    pbool.at(p)=1;
                }
            }

            Est es_(mo,mo1,po,pbool); print_analysis(pr,es,es_);
        } else {
            print_analysis(pr,es,es);
        }
    }

    cout<<"Done!\n";
}
