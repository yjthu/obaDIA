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


void print_data(const Pre& pr)
{
    {
        ofstream ofs("log2_data.txt");
        if (not ofs) throw runtime_error("can't open log2_data.txt");
        ofs<<"Protein";
        if (pr.op().level()>=2) {
            ofs<<"\tPeptide";
            if (pr.op().level()==3) ofs<<"\tFragment";
        }
        for (int j=0; j<pr.op().ncolumns(); j++) ofs<<'\t'<<pr.pheader().at(j);
        for (unsigned p=0;p<pr.yy().size();p++) for (unsigned q=0;q<pr.yy().at(p).size();q++) for (unsigned f=0;f<pr.yy().at(p).at(q).size();f++) {
            ofs<<'\n'<<pr.pidvec().at(p);
            if (pr.op().level()==3) ofs<<'\t'<<pr.qidvec().at(p).at(q);
            if (pr.op().level()>=2) ofs<<'\t'<<pr.fidvec().at(p).at(q).at(f);
            const vector<double>& ref=pr.pqfid_m().find(pr.pidvec().at(p)+'\t'+pr.qidvec().at(p).at(q)+'\t'+pr.fidvec().at(p).at(q).at(f))->second->in_;
            for (int l=0;l<pr.op().ncolumns();l++) {
                if (obs(pr.yy().at(p).at(q).at(f).at(l))) {
                    ofs<<'\t'<<log2(ref.at(l));
                } else {
                    ofs<<"\tNA";
                }
            }
        }
        ofs<<'\n';
    }

    if (pr.op().level()==3) {
        ofstream ofs("fragments_for_protein_quantification.txt");
        if (not ofs) throw runtime_error("can't open fragments_for_protein_quantification.txt");
        ofs<<"Protein";
        if (pr.op().level()>=2) {
            ofs<<"\tPeptide";
            if (pr.op().level()==3) ofs<<"\tFragment";
        }
        for (int j=0; j<pr.op().ncolumns(); j++) ofs<<'\t'<<pr.pheader().at(j);
        for (unsigned p=0;p<pr.yy().size();p++) for (unsigned q=0;q<pr.yy().at(p).size();q++) for (unsigned f=0;f<pr.yy().at(p).at(q).size();f++) {
            ofs<<'\n'<<pr.pidvec().at(p);
            if (pr.op().level()==3) ofs<<'\t'<<pr.qidvec().at(p).at(q);
            if (pr.op().level()>=2) ofs<<'\t'<<pr.fidvec().at(p).at(q).at(f);
            const vector<double>& ref=pr.pqfid_m().find(pr.pidvec().at(p)+'\t'+pr.qidvec().at(p).at(q)+'\t'+pr.fidvec().at(p).at(q).at(f))->second->in_;
            for (int l=0;l<pr.op().ncolumns();l++) {
                if (obs(pr.yy().at(p).at(q).at(f).at(l))) {
                    ofs<<'\t'<<ref.at(l);
                } else {
                    ofs<<"\tNA";
                }
            }
        }
        for (map<string,Row*>::const_iterator it=pr.pqfid_m0().begin();it!=pr.pqfid_m0().end();it++) if (pr.pqfid_m().find(it->first)==pr.pqfid_m().end()) {
            ofs<<'\n'<<it->first;
            for (int l=0;l<pr.op().ncolumns();l++) ofs<<"\tNA";
        }
        ofs<<'\n';
    }


    if (pr.op().level()==3) {
        ofstream ofs("peptide_level.txt");
        if (not ofs) throw runtime_error("can't open peptide_level.txt");
        ofs<<"Protein\tPeptide";
        for (int j=0; j<pr.op().ncolumns(); j++) ofs<<'\t'<<pr.pheader().at(j);
        ofs<<"\tnFragment";
        for (unsigned p=0;p<pr.yy().size();p++) for (unsigned q=0;q<pr.yy().at(p).size();q++) {
            vector<double> qsum(pr.op().ncolumns());
            for (unsigned f=0;f<pr.yy().at(p).at(q).size();f++) {
                const vector<double>& ref=pr.pqfid_m().find(pr.pidvec().at(p)+'\t'+pr.qidvec().at(p).at(q)+'\t'+pr.fidvec().at(p).at(q).at(f))->second->in_;
                for (int l=0;l<pr.op().ncolumns();l++) {
                    qsum.at(l) += (obs(pr.yy().at(p).at(q).at(f).at(l)) ? ref.at(l) : 0);
                }
            }
            ofs<<'\n'<<pr.pidvec().at(p)<<'\t'<<pr.qidvec().at(p).at(q);
            for (int l=0;l<pr.op().ncolumns();l++) ofs<<'\t'<<qsum.at(l);
            ofs<<'\t'<<pr.yy().at(p).at(q).size();
        }
        ofs<<'\n';
    }

    if (pr.op().level()>=2) {
        ofstream ofs("protein_level.txt");
        if (not ofs) throw runtime_error("can't open protein_level.txt");
        ofs<<"Protein";
        for (int j=0; j<pr.op().ncolumns(); j++) ofs<<'\t'<<pr.pheader().at(j);
        if (pr.op().level()==3) ofs<<"\tnFragment";
        ofs<<"\tnPeptide";
        for (unsigned p=0;p<pr.yy().size();p++) {
            multimap<double,int,greater<double> > qm_cor;
            for (unsigned q=0;q<pr.yy().at(p).size();q++) {
                double tmp = accumulate(pr.ym_cor().at(p).at(q).begin(),pr.ym_cor().at(p).at(q).end(),0.);
                tmp /= pr.ym_cor().at(p).at(q).size();
                qm_cor.insert(make_pair(tmp,q));
            }
            multimap<double,int>::const_iterator it=qm_cor.begin();
            vector<int> qindex;
            for (unsigned q=0;static_cast<int>(q)<pr.op().top_q() and q<pr.yy().at(p).size();q++,it++) {
                qindex.push_back(it->second);
            }
            vector<double> psum(pr.op().ncolumns());
            int nfrag=0;
            for (vector<int>::const_iterator q=qindex.begin();q!=qindex.end();q++) for (unsigned f=0;f<pr.yy().at(p).at(*q).size();f++) {
                const vector<double>& ref=pr.pqfid_m().find(pr.pidvec().at(p)+'\t'+pr.qidvec().at(p).at(*q)+'\t'+pr.fidvec().at(p).at(*q).at(f))->second->in_;
                for (int l=0;l<pr.op().ncolumns();l++) {
                    psum.at(l) += (obs(pr.yy().at(p).at(*q).at(f).at(l)) ? ref.at(l) : 0);
                }
                nfrag++;
            }
            ofs<<'\n'<<pr.pidvec().at(p);
            for (int l=0;l<pr.op().ncolumns();l++) ofs<<'\t'<<psum.at(l);
            ofs<<'\t'<<nfrag;
            if (pr.op().level()==3) ofs<<'\t'<<qindex.size();
        }
        ofs<<'\n';
    }




}


void print_selection(const Pre& pr)
{
    ofstream ofs("fragment_selection.txt");
    if (not ofs) throw runtime_error("can't open fragment_selection.txt");
    ofs<<"Protein";
    if (pr.op().level()>=2) {
        ofs<<"\tPeptide";
        if (pr.op().level()==3) ofs<<"\tFragment";
        if (pr.op().shared_q()) ofs<<"\tSHARED_PEPTIDE";
        ofs<<"\tSDF";
    }
    for (int t=0;t<pr.op().nt();t++) ofs<<"\tMIN_OBS_"<<t+1;
    if (pr.op().level()>=2) {
        ofs<<"\tMIN_CORREL";
        if (pr.op().level()==2) ofs<<"\tMIN_PEP_PER_PROT\tMAX_PEP_PER_PROT";
        else ofs<<"\tMIN_FRAG_PER_PEP\tMAX_FRAG_PER_PEP\tMIN_PEP_PER_PROT";
    }
    for (list<Row>::const_iterator it=pr.pdata0().begin();it!=pr.pdata0().end();it++) {
        ofs<<'\n'<<it->pid;
        if (pr.op().level()==3) ofs<<'\t'<<it->qid;
        if (pr.op().level()>=2) {
            ofs<<'\t'<<it->fid;
            if (pr.op().shared_q()) ofs<<'\t'<<(it->pcount>1 ? 'y':'n');
            ofs<<'\t'<<(it->outlier ? 'y':'n');
        }
        for (int t=0;t<pr.op().nt();t++) ofs<<'\t'<<(it->nobs.at(t)<pr.op().min_obs().at(t) ? 'y':'n');
        if (pr.op().level()>=2) {
            ofs<<'\t'<<(it->m_cor<pr.op().min_correl() ? 'y':'n');
            ofs<<'\t'<<(it->fcount<pr.op().min_f() ? 'y':'n')
                <<'\t'<<(it->mcrank>pr.op().max_f()? 'y':'n');
            if (pr.op().level()==3) ofs<<'\t'<<(it->qcount<pr.op().min_q() ? 'y':'n');
        }
    }
    ofs<<'\n';
}






void print_analysis(const Pre& pr,const Est& es,const Est& es_)
{
    const Option& op(pr.op());
    const Post& po(es.po());

    ofstream ofs("analysis_output.txt");
    if (not ofs) throw runtime_error("can't open analysis_output.txt");
    ofs<<"Protein";
    if (op.level()>=2) {
        ofs<<"\tnPeptide";
        if (op.level()==3) ofs<<"\tnFragment";
    }
    ofs<<"\tLabel\tLabel2\tlog2FC\tlog2FC_SE\tscore\tSignedScore\tFDR\tlog_oddsDE";
    if (es_.mo().modulebool()) ofs<<"\tscore_\tFDR_\tlog_oddsDE_";

    if (op.rep()) {
        for (int s=0;s<op.ssize().at(0);s++) ofs<<"\tlog2FC_"<<s+1;
        ofs<<"\tnUp\tnDown";
    }
    vector<vector<int> > cp(op.nc(),vector<int>(po.yy().size(),-1));
    for (int c=0,i=0,j=0;c<op.nc();c++) {
        const int t1=op.ct1().at(c),t2=op.ct2().at(c);
        for (unsigned p=0;p<po.yy().size();p++) if (po.pq().at(c).at(p)==1) {
            ofs<<'\n'<<pr.pidvec().at(p);
            if (op.level()==3) ofs<<'\t'<<count_if(po.nf().at(c).at(p).begin(),po.nf().at(c).at(p).end(),bind2nd(greater_equal<int>(),op.min_f()));
            if (op.level()>=2) ofs<<'\t'<<accumulate(po.nf().at(c).at(p).begin(),po.nf().at(c).at(p).end(),0);
            const double score=1/(exp(-es.lo().at(i))+1);
            ofs<<'\t'<<t1<<'/'<<t2
                <<'\t'<<op.labels().at(t1)<<'/'<<op.labels().at(t2)
                <<'\t'<<po.log2fc().at(c).at(p)
                <<'\t'<<po.log2fc_SE().at(c).at(p)
                <<'\t'<<score
                <<'\t'<<(po.log2fc().at(c).at(p)>0?1:-1)*score
                <<'\t'<<es.fdr().find(es.lo().at(i))->second.fdr<<'\t'<<es.lo().at(i);
            cp.at(c).at(p)=i;
            i++;
            if (es_.mo().modulebool()) {
                if (es_.pbool().at(p)==1) {
                    ofs<<'\t'<<1/(exp(-es_.lo().at(j))+1)<<'\t'<<es_.fdr().find(es_.lo().at(j))->second.fdr<<'\t'<<es_.lo().at(j);
                    j++;
                } else {
                    ofs<<"\tNA\tNA\tNA";
                }
            }
            if (op.rep()) {
                const vector<double>& ref=po.log2fc_s().at(c).at(p);
                for (int s=0;s<op.ssize().at(0);s++) {
                    if (obs(ref.at(s))) ofs<<'\t'<<ref.at(s); else ofs<<"\tNA";
                }
                vector<double> vlog2fc;
                for (unsigned s=0;s<ref.size();s++) {
                    if (obs(ref.at(s))) vlog2fc.push_back(ref.at(s));
                }
                ofs<<'\t'<<count_if(vlog2fc.begin(),vlog2fc.end(),bind2nd(greater<double>(),0))
                    <<'\t'<<count_if(vlog2fc.begin(),vlog2fc.end(),bind2nd(less<double>(),0));
            }
        }
    }
    ofs<<'\n';

    ofstream ofs_("analysis_output_wide_format.txt");
    if (not ofs_) throw runtime_error("can't open analysis_output_wide_format.txt");
    ofs_<<"Comparison";
    for (int c=0;c<op.nc();c++) {
        const int t1=op.ct1().at(c),t2=op.ct2().at(c);
        for (int i=0;i<4;i++) ofs_<<'\t'<<op.labels().at(t1)<<'/'<<op.labels().at(t2);
    }
    ofs_<<"\nProtein";
    for (int c=0;c<op.nc();c++) ofs_<<"\tlog2FC\tscore\tFDR\tlog_oddsDE";

    for (unsigned p=0;p<po.yy().size();p++) {
        ofs_<<'\n'<<pr.pidvec().at(p);
        for (int c=0;c<op.nc();c++) {
            if (po.pq().at(c).at(p)==1) {
                ofs_<<'\t'<<po.log2fc().at(c).at(p)
                    <<'\t'<<1/(exp(-es.lo().at(cp.at(c).at(p)))+1)<<'\t'<<es.fdr().find(es.lo().at(cp.at(c).at(p)))->second.fdr<<'\t'<<es.lo().at(cp.at(c).at(p));
            } else {
                ofs_<<"\t\t\t\t";
            }
        }
    }
    ofs_<<'\n';
}







