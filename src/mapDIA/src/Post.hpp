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
class Pre;


class Post
{
    const Option d_op;
    vector<vector<vector<vector<double> > > > d_zz;
    const vector<vector<vector<vector<double> > > > d_yy;
    const vector<vector<vector<vector<int> > > > d_im;
    vector<vector<vector<double> > > d_log2fc_s;
    vector<vector<double> > d_logf0,d_logf1,d_log2fc,d_log2fc_SE;
    vector<vector<vector<int> > > d_nf;
    vector<vector<int> > d_pq;
    const vector<vector<vector<vector<int> > > > d_yobs;

    void tmp_remove_pair(const int t1,const int t2);
    void tmp_remove(const int t1,const int t2);
    void center_c(const int t1,const int t2);
    void set_nf(const int t1,const int t2,const int c);
    void sigE_c(const int t1,const int t2,double& aE,double& bE);
    void logf_c(const int t1,const int t2,const int c,const double aE,const double bE);
    void log2fc_c(const int t1,const int t2,const int c);


    public:
    Post(const Pre& pr);
    const vector<vector<vector<vector<double> > > >& yy() const { return d_yy; }
    const vector<vector<vector<vector<int> > > >& im() const { return d_im; }
    const vector<vector<vector<vector<double> > > >& zz() const { return d_zz; }
    const vector<vector<double> >& log2fc() const { return d_log2fc; }
    const vector<vector<vector<double> > >& log2fc_s() const { return d_log2fc_s; }
    const vector<vector<double> >& log2fc_SE() const { return d_log2fc_SE; }
    const Option& op() const { return d_op; }
    const vector<vector<double> >& logf0() const { return d_logf0; }
    const vector<vector<double> >& logf1() const { return d_logf1; }
    const vector<vector<int> >& pq() const { return d_pq; }
    const vector<vector<vector<int> > >& nf() const { return d_nf; }
};
