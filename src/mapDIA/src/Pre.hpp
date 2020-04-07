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


class Option;


struct Row
{
    string pid, qid, fid;//protein id, peptide id, fragment id
    vector<double> in;//intensity //pdata0-raw intensity
    vector<double> in_;//raw intensity,pdata only
    vector<int> im;//impute status
    vector<int> nobs;
    double pseudoCV;
    double m_cor;
    double rt;//retention time
    int fcount;//number of fragments from the same peptide
    int qcount;//number of peptides from the same protein
    int pcount;//number of proteins this peptide belong to
    bool outlier;
    int mcrank;
    Row(const int ncolumns,const vector<int>& freq_cut,const double CV_cut,const double mc_cut,const int min_f,const int q_cut) : in(vector<double>(ncolumns)), nobs(freq_cut), pseudoCV(CV_cut), m_cor(mc_cut), rt(0), fcount(min_f), qcount(q_cut), pcount(0), outlier(false), mcrank(0) { }
};


class Pre
{
    map<string,Row*> d_pqfid_m0;
    map<string,Row*> d_pqfid_m;
    list<Row> d_pdata,d_pdata0;
    vector<string> d_pheader,d_pidvec,d_pidvecrm;
    set<string> d_pidset,d_pidset0;
    vector<vector<vector<double> > > d_ym_cor;
    vector<vector<vector<string> > > d_fidvec;
    vector<vector<string> > d_qidvec;
    vector<vector<vector<vector<int> > > > d_yobs;
    vector<vector<vector<vector<double> > > > d_yy;
    vector<vector<vector<vector<int> > > > d_im;
    const Option d_op;

    void round_RT();
    void center_frag();
    void setOutlier();
    void setPseudoCV();
    void setMc();
    void setFcount();
    void setQcount();
    void setShared_q();
    void norm_RT();
    void norm_TIS();
    void filter();
    void setAll();
    void dup();
    void read();
    void preprocess();


    public:
    Pre(const Option& op);
    const map<string,Row*>& pqfid_m0() const { return d_pqfid_m0; }
    const map<string,Row*>& pqfid_m() const { return d_pqfid_m; }
    const set<string>& pidset() const { return d_pidset; }
    const list<Row>& pdata() const { return d_pdata; }
    const list<Row>& pdata0() const { return d_pdata0; }
    const vector<vector<vector<double> > >& ym_cor() const { return d_ym_cor; }
    const vector<vector<vector<string> > >& fidvec() const { return d_fidvec; }
    const vector<vector<string> >& qidvec() const { return d_qidvec; }
    const vector<string>& pheader() const { return d_pheader; }
    const vector<string>& pidvecrm() const { return d_pidvecrm; }
    const vector<string>& pidvec() const { return d_pidvec; }
    const set<string>& pidset0() const { return d_pidset0; }
    const vector<vector<vector<vector<double> > > >& yy() const { return d_yy; }
    const vector<vector<vector<vector<int> > > >& im() const { return d_im; }
    const vector<vector<vector<vector<int> > > >& yobs() const { return d_yobs; }
    const Option& op() const { return d_op; }
};
