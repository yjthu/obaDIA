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


class Module;
class Post;

struct fdrs { double fdr; int t; };


class Est
{
    const Module d_mo,d_mo1;
    const Post d_po;
    const vector<int> d_pbool;
    vector<vector<vector<double> > > d_logp;
    vector<vector<int> > d_Z;
    vector<double> d_Phi,d_lo;//log_odds
    map<double,fdrs,greater<double> > d_fdr;

    int ICM();
    void setFdr();
    void Construct();


    public:
    Est(const Post& po);
    Est(const Module& mo,const Module& mo1,const Post& po,const vector<int>& pbool);
    const Post& po() const { return d_po; }
    const Module& mo() const { return d_mo; }
    const Module& mo1() const { return d_mo1; }
    const vector<vector<vector<double> > >& logp() const { return d_logp; }
    const vector<vector<int> >& Z() const { return d_Z; }
    const vector<double>& Phi() const { return d_Phi; }
    const vector<double>& lo() const { return d_lo; }
    const vector<int>& pbool() const { return d_pbool; }
    const map<double,fdrs,greater<double> >& fdr() const { return d_fdr; }
};
