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


typedef multimap<string,Row*>::const_iterator mmit;
typedef multimap<string,const Row*>::const_iterator cmmit;
const bool rtplot=0;



double median(const vector<double>& data)
{
    vector<double> sorted_data(data);
    sort(sorted_data.begin(),sorted_data.end());
    const size_t lhs = (data.size() - 1) / 2 ;
    const size_t rhs = data.size() / 2 ;
    if (data.size() == 0) throw runtime_error("median");
    if (lhs == rhs) return sorted_data.at(lhs) ;
    else return (sorted_data.at(lhs) + sorted_data.at(rhs))/2.0;
}


class pred_shared_q
{
    const set<string> d_inc;
    public:
    pred_shared_q(const set<string>& inc) : d_inc(inc) { }
    bool operator() (const Row& tmprow) const { return tmprow.pcount>1 and d_inc.find(tmprow.pid)==d_inc.end(); }
};


class pred_freq
{
    const int d_nt;
    const vector<int> d_freq_cut;
    const set<string> d_inc;
    public:
    pred_freq(const int nt,const vector<int>& freq_cut,const set<string>& inc) : d_nt(nt),d_freq_cut(freq_cut),d_inc(inc) { }
    bool operator() (const Row& tmprow) const {
        int count=0;
        for (int t=0;t<d_nt;t++) if (tmprow.nobs.at(t)>=d_freq_cut.at(t)) count++;
        return count<2 and d_inc.find(tmprow.pid)==d_inc.end();
    }
};


class pred_mc
{
    const double d_mc_cut,d_cv_cut;
    const set<string> d_inc;
    public:
    pred_mc(const double mc_cut,const double cv_cut,const set<string>& inc) : d_mc_cut(mc_cut),d_cv_cut(cv_cut),d_inc(inc) { }
    bool operator() (const Row& tmprow) const { return tmprow.pseudoCV>=d_cv_cut and tmprow.m_cor<d_mc_cut and d_inc.find(tmprow.pid)==d_inc.end(); }
    //bool operator() (const Row& tmprow) const { return tmprow.m_cor<d_mc_cut and d_inc.find(tmprow.pid)==d_inc.end(); }
};


class pred_fmin
{
    const int d_min_f;
    const set<string> d_inc;
    public:
    pred_fmin(const int min_f,const set<string>& inc) : d_min_f(min_f),d_inc(inc) { }
    bool operator() (const Row& tmprow) const { return tmprow.fcount<d_min_f and d_inc.find(tmprow.pid)==d_inc.end(); }
};


class pred_fmax
{
    const int d_max_f;
    const set<string> d_inc;
    public:
    pred_fmax(const int max_f,const set<string>& inc) : d_max_f(max_f),d_inc(inc) { }
    bool operator() (const Row& tmprow) const { return tmprow.mcrank>d_max_f and d_inc.find(tmprow.pid)==d_inc.end(); }
};


class pred_qcount
{
    const int d_q_cut;
    const set<string> d_inc;
    public:
    pred_qcount(const int q_cut,const set<string>& inc) : d_q_cut(q_cut),d_inc(inc) { }
    bool operator() (const Row& tmprow) const { return tmprow.qcount<d_q_cut and d_inc.find(tmprow.pid)==d_inc.end(); }
};


void Pre::round_RT()
{
    const double sc = pow(10.,op().dp());
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        it->rt *= sc;
        it->rt = round(it->rt);
        it->rt /= sc;
    }
}


void Pre::center_frag()
{
    if (op().rep()) {
        for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) for (int s=0;s<op().ssize().at(0);s++) {
            vector<double> tmpvec;
            for (int t=0;t<op().nt();t++) if (obs(it->in.at(op().np().at(t)+s))) {
                tmpvec.push_back(it->in.at(op().np().at(t)+s));
            }
            if (tmpvec.size()>1) {
                const double center = median(tmpvec);
                for (int t=0;t<op().nt();t++) if (obs(it->in.at(op().np().at(t)+s))) {
                    it->in.at(op().np().at(t)+s) -= center;
                }
            } else {
                for (int t=0;t<op().nt();t++) it->in.at(op().np().at(t)+s) = NAN;
            }
        }
    } else {
        for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
            vector<double> tmpvec;
            for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) {
                tmpvec.push_back(it->in.at(l));
            }
            if (tmpvec.size()>1) {
                const double center = median(tmpvec);
                for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) {
                    it->in.at(l) -= center;
                }
            } else {
                for (int l=0;l<op().ncolumns();l++) it->in.at(l) = NAN;
            }
        }
    }
}


void Pre::setOutlier()
{
    multimap<string,Row*> pid_mm;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        pid_mm.insert(make_pair(it->pid,&(*it)));
    }
    for (mmit end1,itp=pid_mm.begin();itp!=pid_mm.end();itp=end1) {
        end1=pid_mm.upper_bound(itp->first);
        vector<vector<double> > ww;
        for (mmit it=itp;it!=end1;it++) ww.push_back(it->second->in);
        if (ww.size()==1) continue;
        vector<double> mtrace(op().ncolumns());
        double var=0;
        int nobs=0;
        for (int l=0;l<op().ncolumns();l++) {
            vector<double> tmpvec;
            for (unsigned f=0;f<ww.size();f++) if (obs(ww.at(f).at(l))) {
                tmpvec.push_back(ww.at(f).at(l));
            }
            if (tmpvec.size()>0) mtrace.at(l)=median(tmpvec); else mtrace.at(l)=NAN;
            if (tmpvec.size()>1) {
                const double wbar=accumulate(tmpvec.begin(),tmpvec.end(),0.)/tmpvec.size();
                for (unsigned i=0;i<tmpvec.size();i++) var += pow(tmpvec.at(i)-wbar,2);
                nobs+=tmpvec.size();
            }
        }
        if (nobs>1) {
            var /= nobs-1;
            for (int l=0;l<op().ncolumns();l++) if (obs(mtrace.at(l))) {
                const double tmpub = mtrace.at(l) + op().SDF()*sqrt(var);
                const double tmplb = mtrace.at(l) - op().SDF()*sqrt(var);
                for (mmit it=itp;it!=end1;it++) {
                    double& ref = it->second->in.at(l);
                    if (obs(ref) and (ref>tmpub or ref<tmplb)) {
                        ref=NAN;
                        it->second->outlier=true;
                    }
                }
            }
        }
    }
}


void Pre::setPseudoCV()
{
    multimap<string,Row*> pid_mm;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        pid_mm.insert(make_pair(it->pid,&(*it)));
    }
    for (mmit end1,itp=pid_mm.begin();itp!=pid_mm.end();itp=end1) {
        end1=pid_mm.upper_bound(itp->first);
        double aveCV=0;
        int nf=0;
        for (mmit it=itp;it!=end1;it++) {
            vector<double> in_f;
            for (int l=0;l<op().ncolumns();l++) if (obs(it->second->in.at(l))) {
                in_f.push_back(it->second->in_.at(l));
            }
            const double ave=accumulate(in_f.begin(),in_f.end(),0.)/in_f.size();
            double var=0;
            for (unsigned l=0;l<in_f.size();l++) var+=pow(in_f.at(l)-ave,2)/(in_f.size()-1);
            aveCV+=sqrt(var)/ave;
            nf++;
        }
        aveCV /= nf;
        for (mmit it=itp;it!=end1;it++) it->second->pseudoCV=aveCV;
    }
}


void Pre::setMc()
{
    multimap<string,Row*> pid_mm;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        pid_mm.insert(make_pair(it->pid,&(*it)));
    }
    for (mmit end1,itp=pid_mm.begin();itp!=pid_mm.end();itp=end1) {
        end1=pid_mm.upper_bound(itp->first);
        vector<vector<double> > ww;
        for (mmit it=itp;it!=end1;it++) ww.push_back(it->second->in);
        if (ww.size()==1) {
            itp->second->m_cor=1;
            continue;
        }
        mmit it=itp;
        for (unsigned f=0;f<ww.size();f++,it++) {
            vector<double> cor;
            for (unsigned g=0;g<ww.size();g++) if (f != g) {
                vector<double> center(2);
                vector<double> sd(center.size());
                int ncol0=0;
                for (int l=0;l<op().ncolumns();l++) if (obs(ww.at(f).at(l)) and obs(ww.at(g).at(l))) {
                    center.at(0) += ww.at(f).at(l);
                    center.at(1) += ww.at(g).at(l);
                    ncol0++;
                }
                if (ncol0<2) continue;
                center.at(0) /= ncol0;
                center.at(1) /= ncol0;
                for (int l=0;l<op().ncolumns();l++) if (obs(ww.at(f).at(l)) and obs(ww.at(g).at(l))) {
                    sd.at(0) += pow(ww.at(f).at(l)-center.at(0),2);
                    sd.at(1) += pow(ww.at(g).at(l)-center.at(1),2);
                }
                sd.at(0) = sqrt(sd.at(0));
                sd.at(1) = sqrt(sd.at(1));
                double cortmp=0;
                for (int l=0;l<op().ncolumns();l++) if (obs(ww.at(f).at(l)) and obs(ww.at(g).at(l))) {
                    cortmp += (ww.at(f).at(l)-center.at(0))*(ww.at(g).at(l)-center.at(1));
                }
                if (sd.at(0)==0 or sd.at(1)==0) continue;
                cortmp /= sd.at(0);
                cortmp /= sd.at(1);
                cor.push_back(cortmp);
            }
            if (cor.size()==0) {
                it->second->m_cor=-1;
            } else {
                it->second->m_cor=median(cor);
            }
        }
    }
}


void Pre::setFcount()
{
    multimap<string,Row*> pqid_mm;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        pqid_mm.insert(make_pair(it->pid+'\t'+it->qid,&(*it)));
    }
    for (mmit end1,itq=pqid_mm.begin();itq!=pqid_mm.end();itq=end1) {
        end1=pqid_mm.upper_bound(itq->first);
        int nf=0;
        for (mmit it=itq;it!=end1;it++,nf++);
        for (mmit it=itq;it!=end1;it++) it->second->fcount=nf;

        multimap<double,Row*,greater<double> > mc_mm;
        for (mmit it=itq;it!=end1;it++) {
            mc_mm.insert(make_pair(it->second->m_cor,it->second));
        }
        int m_cor_rank=1;
        for (multimap<double,Row*>::iterator end0,itm=mc_mm.begin();itm!=mc_mm.end();itm=end0) {
            end0=mc_mm.upper_bound(itm->first);
            multimap<double,int*,greater<double> > in_mm;
            for (multimap<double,Row*>::iterator it=itm;it!=end0;it++) {
                double tmp=0;
                for (int l=0;l<op().ncolumns();l++) {
                    if (obs(it->second->in.at(l))) tmp += (it->second->in_.at(l));
                }
                in_mm.insert(make_pair(tmp,&it->second->mcrank));
            }
            for (multimap<double,int*>::iterator it=in_mm.begin();it!=in_mm.end();it++) {
                *it->second=m_cor_rank++;
            }
        }


    }
}


void Pre::setQcount()
{
    multimap<string,Row*> pid_mm;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        pid_mm.insert(make_pair(it->pid,&(*it)));
    }
    for (mmit end1,itp=pid_mm.begin();itp!=pid_mm.end();itp=end1) {
        end1=pid_mm.upper_bound(itp->first);
        set<string> qidset;
        for (mmit it=itp;it!=end1;it++) qidset.insert(it->second->qid);
        for (mmit it=itp;it!=end1;it++) it->second->qcount=qidset.size();
    }
}


void Pre::setShared_q()
{
    multimap<string,Row*> qid_mm;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        if (op().level()==3) {
            qid_mm.insert(make_pair(it->qid,&(*it)));
        } else {
            qid_mm.insert(make_pair(it->fid,&(*it)));
        }
    }
    for (mmit end1,itq=qid_mm.begin();itq!=qid_mm.end();itq=end1) {
        end1=qid_mm.upper_bound(itq->first);
        set<string> pidset;
        for (mmit it=itq;it!=end1;it++) pidset.insert(it->second->pid);
        for (mmit it=itq;it!=end1;it++) it->second->pcount=pidset.size();
    }
}


void Pre::norm_RT()
{
    cout<<"RT("<<op().window()<<") normalization...\n"<<flush;
    vector<double> rsum(pdata().size());
    int i=0;
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++,i++) {
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) rsum.at(i)+=it->in.at(l);
    }
    map<double,vector<double> > denom;
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) denom[it->rt].assign(op().ncolumns(),0);
    map<double,vector<double> > rtsum(denom);
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
        vector<double>& ref=rtsum.find(it->rt)->second;
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) {
            ref.at(l) += it->in.at(l);
        }
    }
    i=0;
    for (map<double,vector<double> >::iterator itd=denom.begin();itd!=denom.end();itd++,i++) {
        if (i%1000==0) cout<<'\r'<<i*100/denom.size()<<'%'<<flush;
        for (map<double,vector<double> >::const_iterator itr=rtsum.begin();itr!=rtsum.end();itr++) {
            double exponent = itr->first - itd->first;
            exponent *= -exponent;
            exponent /= 2 * op().window() * op().window();
            const double w=exp(exponent);
            for (int l=0;l<op().ncolumns();l++) itd->second.at(l) += w * (itr->second.at(l));
        }
    }
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        const vector<double>& ref=denom.find(it->rt)->second;
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) it->in.at(l) /= ref.at(l);
    }
    i=0;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++,i++) {
        double rsum_=0;
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) rsum_+=it->in.at(l);
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) it->in.at(l)*=rsum.at(i)/rsum_;
    }
    cout<<"\r100%\n";
}


void Pre::norm_TIS()
{
    cout<<"TIS normalization... "<<flush;
    double sum=0;
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) sum+=it->in.at(l);
    }
    vector<double> colSums(op().ncolumns());
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) {
        colSums.at(l) += it->in.at(l);
    }
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) {
        it->in.at(l) /= colSums.at(l);
    }
    double sum_norm=0;
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) sum_norm+=it->in.at(l);
    }
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) it->in.at(l)*=sum/sum_norm;
    }
    cout<<"done\n";
}


void Pre::read()
{
    string str0;
    ifstream prot_ifs(op().file().c_str());
    if (not prot_ifs) throw runtime_error("can't open "+op().file());
    getline(prot_ifs,str0);
    istringstream iss(str0);
    getline(iss,str0,'\t');
    if (op().level()>=2) getline(iss,str0,'\t');
    if (op().level()==3) getline(iss,str0,'\t');
    for (int j=0; j<op().ncolumns(); j++) {
        if (not getline(iss,str0,'\t')) throw runtime_error("unequal columns");
        d_pheader.push_back(str0);
    }
    int mRT=0;
    for (int i=2;getline(prot_ifs,str0) and (not str0.empty());i++) {
        Row tmp_f(op().ncolumns(),op().min_obs(),op().min_CV(),op().min_correl(),op().min_f(),op().min_q());
        iss.clear(); iss.str(str0+"\t");
        getline(iss,tmp_f.pid,'\t');
        if (op().level()==2) { tmp_f.qid=tmp_f.pid; getline(iss,tmp_f.fid,'\t'); }
        if (op().level()==3) { getline(iss,tmp_f.qid,'\t'); getline(iss,tmp_f.fid,'\t'); }
        for (int j=0;j<op().ncolumns();j++) {
            if (not getline(iss,str0,'\t')) throw runtime_error("unequal columns: row "+to<string>(i));
            if (str0=="NA" or str0=="0" or str0.empty()) {
                tmp_f.in.at(j)=NAN;
            } else {
                try {
                    tmp_f.in.at(j)=(op().log2()?to<double>(str0):pow(2.,to<double>(str0)));
                } catch (runtime_error& e) {
                    cerr<<e.what()<<": row "<<i<<", column "<<j+1+op().level();
                    exit(1);
                }
            }
        }
        if (rtplot or op().normalization()=="RT") {
            if (not getline(iss,str0,'\t')) throw runtime_error("RT time: row "+to<string>(i));
            if (str0=="NA" or str0.empty()) {
                mRT++; 
            } else {
                try {
                    tmp_f.rt=to<double>(str0);
                } catch (runtime_error& e) {
                    cerr<<e.what()<<" RT time : row "<<i;
                    exit(1);
                }
            }
        }
        d_pdata.push_back(tmp_f);
    }
    if (mRT>0) cout<<mRT<<" Row(s) with missing retention time data are removed\n";
}


void Pre::preprocess()
{
    dup();

    if (op().normalization()=="RT") {
        if (op().dp()>INT_MIN) round_RT(); //round RT to shorten RT normalisation
        norm_RT();
    } else if (op().normalization()=="TIS") {
        norm_TIS();
    } else {
        cout<<"no normalisation\n";
    }


    d_pdata0=pdata();
    for (list<Row>::const_iterator it=pdata0().begin();it!=pdata0().end();it++) {
        d_pidset0.insert(it->pid);
    }
    for (list<Row>::iterator it=d_pdata0.begin();it!=d_pdata0.end();it++) {
        d_pqfid_m0.insert(make_pair(it->pid+'\t'+it->qid+'\t'+it->fid,&(*it)));
    }
    if (pqfid_m0().size()<pdata0().size()) throw runtime_error("duplicated fragment data");

    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        it->in_=it->in;
        for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) it->in.at(l)=log2(it->in.at(l));
    }
}


void Pre::filter()
{
    multimap<string,Row*> pid_mm_;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        pid_mm_.insert(make_pair(it->pid,&(*it)));
    }
    if (op().level()>=2 and op().shared_q()) {
        setShared_q();
        for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
            d_pqfid_m0[it->pid+'\t'+it->qid+'\t'+it->fid]->pcount = it->pcount;
        }
        d_pdata.remove_if(pred_shared_q(op().inclusion()));
    }

    if (op().SDF()>0) {
        setOutlier();
        for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
            d_pqfid_m0[it->pid+'\t'+it->qid+'\t'+it->fid]->outlier = it->outlier;
        }
    }

    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        it->nobs.assign(op().nt(),0);
        for (int t=0;t<op().nt();t++) for (int s=0;s<op().ssize().at(t);s++) if (obs(it->in.at(op().np().at(t)+s))) {
            it->nobs.at(t)++;
        }
    }
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
        d_pqfid_m0[it->pid+'\t'+it->qid+'\t'+it->fid]->nobs = it->nobs;
    }
    if (op().nc()>0) d_pdata.remove_if(pred_freq(op().nt(),op().min_obs(),op().inclusion()));

    if (op().level()>=2) {
        setPseudoCV();
        setMc();
        //int ia=0;
        for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
            d_pqfid_m0[it->pid+'\t'+it->qid+'\t'+it->fid]->pseudoCV = it->pseudoCV;
            //if(ia++>-1)cout<<it->pid+'\t'+it->qid+'\t'+it->fid+' '<<it->pseudoCV<<'\n';
            d_pqfid_m0[it->pid+'\t'+it->qid+'\t'+it->fid]->m_cor = it->m_cor;
        }
        d_pdata.remove_if(pred_mc(op().min_correl(),op().min_CV(),op().inclusion()));
    }

    if (op().level()>=2) {
        setFcount();
        for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
            d_pqfid_m0[it->pid+'\t'+it->qid+'\t'+it->fid]->fcount = it->fcount;
            d_pqfid_m0[it->pid+'\t'+it->qid+'\t'+it->fid]->mcrank = it->mcrank;
        }
        d_pdata.remove_if(pred_fmin(op().min_f(),op().inclusion()));
        d_pdata.remove_if(pred_fmax(op().max_f(),op().inclusion()));
    }

    if (op().level()==3) {
        setQcount();
        for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
            d_pqfid_m0[it->pid+'\t'+it->qid+'\t'+it->fid]->qcount = it->qcount;
        }
        d_pdata.remove_if(pred_qcount(op().min_q(),op().inclusion()));
    }

    if (op().impute_scale()<=0) {
        for (set<string>::const_iterator it0=op().inclusion().begin();it0!=op().inclusion().end();it0++) {
            pair <mmit, mmit> ret = pid_mm_.equal_range(*it0);
            for (mmit it1=ret.first; it1!=ret.second; ++it1) {
                Row* it=it1->second;
                double min_in=DBL_MAX;
                for (int l=0;l<op().ncolumns();l++) {
                    if (obs(it->in.at(l)) and it->in_.at(l)<min_in) min_in=it->in_.at(l);
                }
                it->im.assign(op().ncolumns(),0);
                if (min_in<DBL_MAX) {
                    for (int l=0;l<op().ncolumns();l++) {
                        if (not obs(it->in.at(l))) {
                            it->in_.at(l)=min_in*.5;//op().impute_scale();
                            it->im.at(l)=1;
                        }
                    }
                    it->in=it->in_;
                    for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) it->in.at(l)=log2(it->in.at(l));
                }
            }
        }
    }

    if (op().impute_scale()>0) {
        if (op().impute_type()=="row") for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
            double min_in=DBL_MAX;
            for (int l=0;l<op().ncolumns();l++) {
                if (obs(it->in.at(l)) and it->in_.at(l)<min_in) min_in=it->in_.at(l);
            }
            it->im.assign(op().ncolumns(),0);
            if (min_in<DBL_MAX) {
                for (int l=0;l<op().ncolumns();l++) {
                    if (not obs(it->in.at(l))) {
                        it->in_.at(l)=min_in*op().impute_scale();
                        it->im.at(l)=1;
                    }
                }
                it->in=it->in_;
                for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) it->in.at(l)=log2(it->in.at(l));
            }
        }
        if (op().impute_type()=="group") {
            vector<double> colMin(op().ncolumns(),DBL_MAX);
            for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
                for (int l=0;l<op().ncolumns();l++) {
                    if (obs(it->in.at(l)) and it->in_.at(l)<colMin.at(l)) colMin.at(l)=it->in_.at(l);
                }
            }
            for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
                it->im.assign(op().ncolumns(),0);
                for (int t=0;t<op().nt();t++) {
                    double min_in=DBL_MAX;
                    for (int s=0;s<op().ssize().at(t);s++) if (obs(it->in.at(op().np().at(t)+s)) and it->in_.at(op().np().at(t)+s)<min_in) {
                        min_in=it->in_.at(op().np().at(t)+s);
                    }
                    if (min_in<DBL_MAX) {
                        for (int s=0;s<op().ssize().at(t);s++) {
                            if (not obs(it->in.at(op().np().at(t)+s))) {
                                it->in_.at(op().np().at(t)+s)=min_in*op().impute_scale();
                                it->im.at(op().np().at(t)+s)=1;
                            }
                        }
                    } else {
                        for (int s=0;s<op().ssize().at(t);s++) {
                            if (not obs(it->in.at(op().np().at(t)+s))) {
                                it->in_.at(op().np().at(t)+s)=colMin.at(op().np().at(t)+s)*op().impute_scale();
                                it->im.at(op().np().at(t)+s)=1;
                            }
                        }
                    }
                }
                it->in=it->in_;
                for (int l=0;l<op().ncolumns();l++) if (obs(it->in.at(l))) it->in.at(l)=log2(it->in.at(l));
            }
        }
    }
}


void Pre::setAll()
{
    multimap<string,const Row*> pid_mm;
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
        pid_mm.insert(make_pair(it->pid,&(*it)));
    }
    for (list<Row>::const_iterator it=pdata().begin();it!=pdata().end();it++) {
        d_pidset.insert(it->pid);
    }
    for (set<string>::const_iterator it=pidset0().begin();it!=pidset0().end();it++) {
        if (pidset().end()==pidset().find(*it)) d_pidvecrm.push_back(*it);
    }
    d_pidvec.assign(pidset().begin(),pidset().end()); 
    d_qidvec.resize(pidset().size()); 
    d_yy.resize(pidset().size()); 
    d_im.resize(pidset().size()); 
    d_yobs.resize(pidset().size()); 
    d_fidvec.resize(pidset().size()); 
    d_ym_cor.resize(pidset().size()); 
    unsigned p=0;
    for (cmmit end1,itp=pid_mm.begin();itp!=pid_mm.end();itp=end1,p++) {
        end1=pid_mm.upper_bound(itp->first);
        multimap<string,const Row*> qid_mm;
        set<string> qidset;
        for (cmmit it=itp;it!=end1;it++) {
            qid_mm.insert(make_pair(it->second->qid,it->second));
            qidset.insert(it->second->qid);
        }
        d_qidvec.at(p).assign(qidset.begin(),qidset.end());
        d_yy.at(p).resize(qidvec().at(p).size());
        d_im.at(p).resize(qidvec().at(p).size());
        d_yobs.at(p).resize(qidvec().at(p).size());
        d_fidvec.at(p).resize(qidvec().at(p).size());
        d_ym_cor.at(p).resize(qidvec().at(p).size());
        unsigned q=0;
        for (cmmit end0,itq=qid_mm.begin();itq!=qid_mm.end();itq=end0,q++) {
            end0=qid_mm.upper_bound(itq->first);
            for (cmmit it=itq;it!=end0;it++) {
                d_fidvec.at(p).at(q).push_back(it->second->fid);
                d_yy.at(p).at(q).push_back(it->second->in);
                d_im.at(p).at(q).push_back(it->second->im);
                d_yobs.at(p).at(q).push_back(it->second->nobs);
                d_ym_cor.at(p).at(q).push_back(it->second->m_cor);
            }
        }
    }
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        d_pqfid_m.insert(make_pair(it->pid+'\t'+it->qid+'\t'+it->fid,&(*it)));
    }
}


Pre::Pre(const Option& op) :
    d_op(op)
{
    read();
    preprocess();
    center_frag();
    filter();
    center_frag();
    setAll();
}


void Pre::dup()
{
    multimap<string,Row*> pqfid_mm;
    for (list<Row>::iterator it=d_pdata.begin();it!=d_pdata.end();it++) {
        pqfid_mm.insert(make_pair(it->pid+'\t'+it->qid+'\t'+it->fid,&(*it)));
    }
    vector<string> duprow;
    for (mmit end1,itpqf=pqfid_mm.begin();itpqf!=pqfid_mm.end();itpqf=end1) {
        end1=pqfid_mm.upper_bound(itpqf->first);
        int ndup=0;
        for (mmit it=itpqf;it!=end1;it++,ndup++) {
            if (ndup>0) {
                if (ndup==1) {
                    string str0=it->second->pid;
                    if (op().level()==3) str0 += '\t'+it->second->qid;
                    if (op().level()>=2) str0 += '\t'+it->second->fid;
                    duprow.push_back(str0);
                }
                if (op().level()>=2) {
                    it->second->fid += "_duplicate"+to<string>(ndup);
                } else {
                    it->second->pid += "_duplicate"+to<string>(ndup);
                }
            }
        }
    }
    if (duprow.size()>0) {
        cout<<"duplicate row names: \"duplicates.txt\"\n";
        ofstream ofs("duplicates.txt");
        if (not ofs) throw runtime_error("can't open duplicates.txt");
        copy(duprow.begin(),duprow.end(),ostream_iterator<string>(ofs,"\n"));
    }
}
