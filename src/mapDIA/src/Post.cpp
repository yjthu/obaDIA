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
#include"Pre.hpp"
#include"Post.hpp"


double quantile(const vector<double>& data,const double f)
{
    vector<double> sorted_data(data);
    sort(sorted_data.begin(),sorted_data.end());
    const double index = f * (data.size() - 1) ;
    const size_t lhs = static_cast<size_t>(index) ;
    const double delta = index - lhs ;
    double result;
    if (data.size() == 0) throw runtime_error("quantile");
    if (lhs == data.size() - 1) {
        result = sorted_data.at(lhs) ;
    } else {
        result = (1 - delta) * sorted_data.at(lhs) + delta * sorted_data.at(lhs + 1) ;
    }
    return result ;
}


void Post::tmp_remove(const int t1,const int t2)
{
    for (unsigned p=0;p<yy().size();p++) for (unsigned q=0;q<yy().at(p).size();q++) for (unsigned f=0;f<yy().at(p).at(q).size();f++) {
        if (d_yobs.at(p).at(q).at(f).at(t1)<op().min_obs().at(t1) or d_yobs.at(p).at(q).at(f).at(t2)<op().min_obs().at(t2)) {
            for (int s=0;s<op().ssize().at(t1);s++) d_zz.at(p).at(q).at(f).at(op().np().at(t1)+s)=NAN;
            for (int s=0;s<op().ssize().at(t2);s++) d_zz.at(p).at(q).at(f).at(op().np().at(t2)+s)=NAN;
        }
    }
}


void Post::center_c(const int t1,const int t2)
{
    for (unsigned p=0;p<yy().size();p++) for (unsigned q=0;q<yy().at(p).size();q++) {
        vector<double> tmpvec;
        for (unsigned f=0;f<yy().at(p).at(q).size();f++) {
            const vector<double>& ref=zz().at(p).at(q).at(f);
            for (int s=0;s<op().ssize().at(t1);s++) if (obs(ref.at(op().np().at(t1)+s))) {
                tmpvec.push_back(ref.at(op().np().at(t1)+s));
            }
            for (int s=0;s<op().ssize().at(t2);s++) if (obs(ref.at(op().np().at(t2)+s))) {
                tmpvec.push_back(ref.at(op().np().at(t2)+s));
            }
        }
        const double center = accumulate(tmpvec.begin(),tmpvec.end(),0.)/tmpvec.size();
        for (unsigned f=0;f<yy().at(p).at(q).size();f++) {
            vector<double>& ref=d_zz.at(p).at(q).at(f);
            for (int s=0;s<op().ssize().at(t1);s++) if (obs(ref.at(op().np().at(t1)+s))) {
                ref.at(op().np().at(t1)+s) -= center;
            }
            for (int s=0;s<op().ssize().at(t2);s++) if (obs(ref.at(op().np().at(t2)+s))) {
                ref.at(op().np().at(t2)+s) -= center;
            }
        }
    }
}


void Post::set_nf(const int t1,const int t2,const int c)
{
    for (unsigned p=0;p<yy().size();p++) {
        d_nf.at(c).at(p).assign(yy().at(p).size(),0);
        for (unsigned q=0;q<yy().at(p).size();q++) {
            for (unsigned f=0;f<yy().at(p).at(q).size();f++) {
                bool flag=true;
                for (int s=0;s<op().ssize().at(t1);s++) if (obs(zz().at(p).at(q).at(f).at(op().np().at(t1)+s))) {
                    d_nf.at(c).at(p).at(q)++;
                    flag=false;
                    break;
                }
                if (flag) for (int s=0;s<op().ssize().at(t2);s++) if (obs(zz().at(p).at(q).at(f).at(op().np().at(t2)+s))) {
                    d_nf.at(c).at(p).at(q)++;
                    break;
                }
            }
        }
    }
}


void Post::sigE_c(const int t1,const int t2,double& aE,double& bE)
{
    vector<double> sigE;
    for (unsigned p=0;p<yy().size();p++) for (unsigned q=0;q<yy().at(p).size();q++) {
        double sigEpq=0;
        int nobs=0;
        for (unsigned f=0;f<yy().at(p).at(q).size();f++) {
            const vector<double>& ref=zz().at(p).at(q).at(f);
            for (int s=0;s<op().ssize().at(t1);s++) if (obs(ref.at(op().np().at(t1)+s))) {
                sigEpq += pow(ref.at(op().np().at(t1)+s),2);
                nobs++;
            }
            for (int s=0;s<op().ssize().at(t2);s++) if (obs(ref.at(op().np().at(t2)+s))) {
                sigEpq += pow(ref.at(op().np().at(t2)+s),2);
                nobs++;
            }
        }
        if (nobs-1>0 and sigEpq>0) sigE.push_back(sigEpq/(nobs-1));
    }
    const double mom1 = accumulate(sigE.begin(),sigE.end(),0.)/sigE.size();
    double mom2 = 0;
    for (unsigned i=0; i<sigE.size(); i++) mom2 += sigE.at(i)*sigE.at(i);
    mom2 /= sigE.size();
    aE=(2*mom2-mom1*mom1)/(mom2-mom1*mom1);//MOM estimates
    bE=mom1*mom2/(mom2-mom1*mom1);         //MOM estimates
}


const double V0=1000;
double mlike(const int n,const double sumsq,const double sum,
        const double aE,const double bE,const double ff)
{
    return -.5*log(1+n*V0)-(aE+n/2.)*log(bE+.5*(sumsq-sum*sum/(1./V0+n)+ff));
}


double mlike(const int n1,const double sumsq1,const double sum1,
        const int n2,const double sumsq2,const double sum2,
        const double aE,const double bE,const double ff)
{
    return -.5*(log(1+n1*V0)+log(1+n2*V0))-(aE+(n1+n2)/2.)*log(bE+.5*(sumsq1-sum1*sum1/(1./V0+n1)+sumsq2-sum2*sum2/(1./V0+n2)+ff));
}


void Post::logf_c(const int t1,const int t2,const int c,const double aE,const double bE)
{
    vector<double> sumsqs;
    if (op().fudge()>=0) {
        for (unsigned p=0;p<yy().size();p++) if (pq().at(c).at(p)==1) {
            for (unsigned q=0;q<yy().at(p).size();q++) {
                double sumsqD1=0,sumsqD2=0;
                int n1=0,n2=0;
                for (unsigned f=0;f<yy().at(p).at(q).size();f++) {
                    const vector<double>& ref=zz().at(p).at(q).at(f);
                    for (int s=0;s<op().ssize().at(t1);s++) if (obs(ref.at(op().np().at(t1)+s)) and (im().at(p).at(q).at(f).empty() or im().at(p).at(q).at(f).at(op().np().at(t1)+s)==0)) {
                        n1++;
                        sumsqD1 += pow(ref.at(op().np().at(t1)+s),2);
                    }
                    for (int s=0;s<op().ssize().at(t2);s++) if (obs(ref.at(op().np().at(t2)+s)) and (im().at(p).at(q).at(f).empty() or im().at(p).at(q).at(f).at(op().np().at(t2)+s)==0)) {
                        n2++;
                        sumsqD2 += pow(ref.at(op().np().at(t2)+s),2);
                    }
                }
                if (nf().at(c).at(p).at(q)<op().min_f() or n1<2 or n2<2) continue;
                if (sumsqD1+sumsqD2!=0) sumsqs.push_back(sumsqD1+sumsqD2);
            }
        }
    }
    const double ff=(op().fudge()<0?0:quantile(sumsqs,op().fudge()));

    for (unsigned p=0;p<yy().size();p++) if (pq().at(c).at(p)==1) {
        for (unsigned q=0;q<yy().at(p).size();q++) {
            double sumsqD1=0,sumD1=0,sumsqD2=0;
            int n1=0,n2=0;
            for (unsigned f=0;f<yy().at(p).at(q).size();f++) {
                const vector<double>& ref=zz().at(p).at(q).at(f);
                for (int s=0;s<op().ssize().at(t1);s++) if (obs(ref.at(op().np().at(t1)+s))) {
                    n1++;
                    sumD1 += ref.at(op().np().at(t1)+s);
                    sumsqD1 += pow(ref.at(op().np().at(t1)+s),2);
                }
                for (int s=0;s<op().ssize().at(t2);s++) if (obs(ref.at(op().np().at(t2)+s))) {
                    n2++;
                    sumsqD2 += pow(ref.at(op().np().at(t2)+s),2);
                }
            }
            if (nf().at(c).at(p).at(q)<op().min_f()) continue;
            if (n1<2 or n2<2) { d_nf.at(c).at(p).at(q)=0; continue; }
            d_logf0.at(c).at(p) += mlike(n1+n2,sumsqD1+sumsqD2,0,aE,bE,ff);
            d_logf1.at(c).at(p) += mlike(n1,sumsqD1,sumD1,n2,sumsqD2,-sumD1,aE,bE,ff);
        }
        if (logf0().at(c).at(p)==0 and logf1().at(c).at(p)==0) d_pq.at(c).at(p)=0;
    }
}


void Post::tmp_remove_pair(const int t1,const int t2)
{
    for (unsigned p=0;p<yy().size();p++) for (unsigned q=0;q<yy().at(p).size();q++) for (unsigned f=0;f<yy().at(p).at(q).size();f++) {
        int npair=0;
        vector<double>& ref=d_zz.at(p).at(q).at(f);
        for (int s=0;s<op().ssize().at(0);s++) {
            if (obs(ref.at(op().np().at(t1)+s)) and obs(ref.at(op().np().at(t2)+s))) {
                npair++;
            } else if (obs(ref.at(op().np().at(t1)+s)) != obs(ref.at(op().np().at(t2)+s))) {
                ref.at(op().np().at(t1)+s)=NAN;
                ref.at(op().np().at(t2)+s)=NAN;
            }
        }
        if (npair<op().min_obs().at(0)) {
            for (int s=0;s<op().ssize().at(0);s++) ref.at(op().np().at(t1)+s)=NAN;
            for (int s=0;s<op().ssize().at(0);s++) ref.at(op().np().at(t2)+s)=NAN;
        }
    }
}


void Post::log2fc_c(const int t1,const int t2,const int c) 
{
    for (unsigned p=0;p<zz().size();p++) if (pq().at(c).at(p)==1) {
        if (op().rep()) {
            vector<double> diffvec;
            for (int s=0;s<op().ssize().at(0);s++) {
                vector<double> D;
                for (unsigned q=0;q<zz().at(p).size();q++) for (unsigned f=0;f<zz().at(p).at(q).size();f++) {
                    const vector<double>& ref=zz().at(p).at(q).at(f);
                    if (obs(ref.at(op().np().at(t1)+s))) {
                        D.push_back(ref.at(op().np().at(t1)+s)-ref.at(op().np().at(t2)+s));
                    }
                }
                copy(D.begin(),D.end(),back_inserter(diffvec));
                d_log2fc_s.at(c).at(p).push_back(D.size()>0 ? accumulate(D.begin(),D.end(),0.)/D.size() : NAN);
            }
            d_log2fc.at(c).at(p)=accumulate(diffvec.begin(),diffvec.end(),0.)/diffvec.size();
            double log2fc_var=0;
            for (unsigned s=0;s<diffvec.size();s++) {
                log2fc_var+=pow(diffvec.at(s)-d_log2fc.at(c).at(p),2);
            }
            d_log2fc_SE.at(c).at(p)=sqrt(log2fc_var/(diffvec.size()-2));
        } else {
            double log2_1=0,log2_2=0;
            int nobs1=0,nobs2=0;
            for (unsigned q=0;q<zz().at(p).size();q++) {
                for (unsigned f=0;f<zz().at(p).at(q).size();f++) {
                    const vector<double>& ref=zz().at(p).at(q).at(f);
                    for (int s=0;s<op().ssize().at(t1);s++) if (obs(ref.at(op().np().at(t1)+s))) {
                        log2_1 += ref.at(op().np().at(t1)+s);
                        nobs1++;
                    }
                    for (int s=0;s<op().ssize().at(t2);s++) if (obs(ref.at(op().np().at(t2)+s))) {
                        log2_2 += ref.at(op().np().at(t2)+s);
                        nobs2++;
                    }
                }
            }
            d_log2fc.at(c).at(p)=(nobs1>0 and nobs2>0 ? log2_1/nobs1-log2_2/nobs2 : 0);
            double s1=0,s2=0;
            for (unsigned q=0;q<zz().at(p).size();q++) for (unsigned f=0;f<zz().at(p).at(q).size();f++) {
                const vector<double>& ref=zz().at(p).at(q).at(f);
                for (int s=0;s<op().ssize().at(t1);s++) if (obs(ref.at(op().np().at(t1)+s))) {
                    s1 += pow(ref.at(op().np().at(t1)+s)-log2_1/nobs1,2);
                }
                for (int s=0;s<op().ssize().at(t2);s++) if (obs(ref.at(op().np().at(t2)+s))) {
                    s2 += pow(ref.at(op().np().at(t2)+s)-log2_2/nobs2,2);
                }
            }
            d_log2fc_SE.at(c).at(p)=sqrt((s1+s2)/(nobs1+nobs2-2));
        }
    }
}


Post::Post(const Pre& pr) :
    d_op(pr.op()),d_zz(pr.yy()),d_yy(pr.yy()),d_im(pr.im()),d_yobs(pr.yobs())
{
    d_logf0.assign(op().nc(),vector<double> (yy().size()));
    d_logf1.assign(op().nc(),vector<double> (yy().size()));
    d_log2fc.assign(op().nc(),vector<double> (yy().size()));
    d_log2fc_s.assign(op().nc(),vector<vector<double> > (yy().size()));
    d_log2fc_SE.assign(op().nc(),vector<double> (yy().size()));
    d_nf.assign(op().nc(),vector<vector<int> >(yy().size()));
    d_pq.assign(op().nc(),vector<int>(yy().size()));
    for (int c=0;c<op().nc();c++) {
        d_zz=d_yy;
        const int t1=op().ct1().at(c),t2=op().ct2().at(c);
        if (d_op.rep()) tmp_remove_pair(t1,t2); else tmp_remove(t1,t2);
        set_nf(t1,t2,c);
        for (unsigned p=0;p<yy().size();p++) {
            if (count_if(nf().at(c).at(p).begin(),nf().at(c).at(p).end(),bind2nd(greater_equal<int>(),op().min_f()))>=op().min_q()) d_pq.at(c).at(p)=1;
        }
        center_c(t1,t2);
        double aE=0,bE=0;
        sigE_c(t1,t2,aE,bE);
        cout<<op().labels().at(t1)<<" vs "<<op().labels().at(t2)<<", a="<<aE<<", b="<<bE<<"\n";
        logf_c(t1,t2,c,aE,bE);
        log2fc_c(t1,t2,c);
    }
}
