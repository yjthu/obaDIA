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
#include"Module.hpp"
#include"Post.hpp"
#include"Est.hpp"
#include"nlopt.hpp"


double logcond_fn(const int p,const vector<double>& Phi,const vector<int>& Zij,const Module& mo,const Module& mo1)
{
    double f0=0,f1=0;
    if (mo.modulebool() and mo.adj().at(p).size()>0) {
        for (set<int>::const_iterator it=mo.adj().at(p).begin();it!=mo.adj().at(p).end();it++) {
            if (mo.mrf()=="0_1") {
                f0 += Zij.at(*it);
            } else {
                f0 += 2*Zij.at(*it)-1;
            }
        }
        f0 /= mo.adj().at(p).size();
        f0 *= Phi.at(1);
    }
    if (mo1.modulebool() and mo1.adj().at(p).size()>0) {
        for (set<int>::const_iterator it=mo1.adj().at(p).begin();it!=mo1.adj().at(p).end();it++) {
            if (mo1.mrf()=="0_1") {
                f1 += Zij.at(*it);
            } else {
                f1 += 2*Zij.at(*it)-1;
            }
        }
        f1 /= mo1.adj().at(p).size();
        f1 *= Phi.at(2);
    }
    const double F = Phi.at(0)+f0+f1;
    return Zij.at(p)*F - log(1 + exp(F));
}


class Logpseudol
{
    const vector<vector<int> > d_Z;
    const Post d_po;
    const Module d_mo,d_mo1;
    const vector<int> d_pbool;
    public:
    Logpseudol(const vector<vector<int> >& Z,const Post& po,const Module& mo,const Module& mo1,const vector<int>& pbool) : d_Z(Z),d_po(po),d_mo(mo),d_mo1(mo1),d_pbool(pbool) { }
    double operator()(const vector<double> &d_Phi, vector<double> &grad) const
    {
        double logpseudo=0;
        for (int c=0;c<d_po.op().nc();c++) for (unsigned p=0;p<d_po.yy().size();p++) if (d_pbool.at(p)==1 and d_po.pq().at(c).at(p)==1) {
            logpseudo += logcond_fn(p,d_Phi,d_Z.at(c),d_mo,d_mo1);
        }
        return logpseudo;
    }
    static double wrap(const vector<double> &x, vector<double> &grad, void *data)
    {
        return (*reinterpret_cast<Logpseudol*>(data))(x, grad); 
    }
};


int Est::ICM()
{
    int sumdiff=0;
    for (int c=0;c<d_po.op().nc();c++) {
        vector<int> Zij(Z().at(c));
        for (unsigned p=0;p<d_po.yy().size();p++) if (pbool().at(p)==1 and d_po.pq().at(c).at(p)==1) {
            vector<int> Zij1(Z().at(c)),Zij0(Z().at(c));
            Zij1.at(p) = 1;
            Zij0.at(p) = 0;
            d_logp.at(c).at(p).at(1) = d_po.logf1().at(c).at(p) + logcond_fn(p,Phi(),Zij1,mo(),mo1());
            d_logp.at(c).at(p).at(0) = d_po.logf0().at(c).at(p) + logcond_fn(p,Phi(),Zij0,mo(),mo1());
            d_Z.at(c).at(p)=(logp().at(c).at(p).at(1)>logp().at(c).at(p).at(0) ? 1:0);
        }
        for (unsigned p=0;p<d_po.yy().size();p++) if (pbool().at(p)==1 and d_po.pq().at(c).at(p)==1) sumdiff += (Zij.at(p)!=Z().at(c).at(p));
    }
    return sumdiff;
}


void Est::setFdr()
{

    for (int c=0;c<d_po.op().nc();c++) for (unsigned p=0;p<d_po.yy().size();p++) if (pbool().at(p)==1 and d_po.pq().at(c).at(p)==1) {
        d_lo.push_back(logp().at(c).at(p).at(1)-logp().at(c).at(p).at(0));
    }


    for (unsigned i=0;i<lo().size();i++) d_fdr[lo().at(i)].t++;
    int denom=fdr().begin()->second.t;
    double numer=(1-1/(exp(-fdr().begin()->first)+1))*denom;
    map<double,fdrs>::iterator it=d_fdr.begin();
    for (it++;it!=d_fdr.end();it++) {
        it->second.fdr=numer/denom;
        numer += (1-1/(exp(-it->first)+1))*it->second.t;
        denom += it->second.t;
    }
}


void Est::Construct()
{
    d_Z.assign(d_po.op().nc(),vector<int>(d_po.yy().size()));
    for (int c=0;c<d_po.op().nc();c++) for (unsigned p=0;p<d_po.yy().size();p++) if ((c+p)&1) d_Z.at(c).at(p)=1;
    d_logp.assign(d_po.op().nc(),vector<vector<double> >(d_po.yy().size(),vector<double>(2)));

    ofstream ofs("param.txt");
    if (not ofs) throw runtime_error("can't open param.txt");
    for (int iter=0,sumdiff=1;iter<100 and sumdiff>0;iter++) {
        d_Phi.assign(1,0);
        if (mo().modulebool()) {
            d_Phi.push_back(1);
            if (mo1().modulebool()) d_Phi.push_back(1);
        }
        nlopt::opt opt(nlopt::LN_COBYLA,Phi().size());
        if (mo().modulebool()) {
            vector<double> lb(Phi().size(),1e-9); lb.at(0) = -10;
            opt.set_lower_bounds(lb);
            opt.set_upper_bounds(vector<double>(Phi().size(),10));
        } else {
            opt.set_lower_bounds(log(d_po.op().minDE()/(1-d_po.op().minDE())));
            opt.set_upper_bounds(log(d_po.op().maxDE()/(1-d_po.op().maxDE())));
        }
        opt.set_maxeval(10000);
        opt.set_xtol_abs(1e-6);
        Logpseudol obj_fn(Z(),d_po,mo(),mo1(),pbool());
        opt.set_max_objective(Logpseudol::wrap, &obj_fn);
        double maxf;
        nlopt::result result = opt.optimize(d_Phi, maxf);
        if (result<0) throw runtime_error("result");

        sumdiff=ICM();

        if (mo1().modulebool()) {
            ofs<<"iteration "<<iter<<": gamma = "<<Phi().at(0)<<" beta1 = "<<Phi().at(1)<<" beta2 = "<<Phi().at(2)<<'\n';
            cout<<"iteration "<<iter<<": gamma = "<<Phi().at(0)<<" beta1 = "<<Phi().at(1)<<" beta2 = "<<Phi().at(2)<<'\n';
        } else if (mo().modulebool()) {
            ofs<<"iteration "<<iter<<": gamma = "<<Phi().at(0)<<" beta = "<<Phi().at(1)<<'\n';
            cout<<"iteration "<<iter<<": gamma = "<<Phi().at(0)<<" beta = "<<Phi().at(1)<<'\n';
        } else {
            ofs<<"iteration "<<iter<<": gamma = "<<Phi().at(0)<<'\t'
                <<"proportion DE = "<<1/(1+exp(-Phi().at(0)))<<'\n';
            cout<<"iteration "<<iter<<": gamma = "<<Phi().at(0)<<'\t'
                <<"proportion DE = "<<1/(1+exp(-Phi().at(0)))<<'\n';
        }
    }

    setFdr();
}


Est::Est(const Post& po) :
    d_mo(Module()),d_mo1(Module()),d_po(po),d_pbool(vector<int> (po.yy().size(),1))
{
    Construct();
}


Est::Est(const Module& mo,const Module& mo1,const Post& po,const vector<int>& pbool) :
    d_mo(mo),d_mo1(mo1),d_po(po),d_pbool(pbool)
{
    Construct();
}
