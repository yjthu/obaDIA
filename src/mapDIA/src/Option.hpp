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


class Option
{
    const map<string,string> d_opm;
    bool d_rep,d_log2,d_shared_q;
    string d_normalization,d_file,d_impute_type;
    int d_min_f,d_max_f,d_min_q,d_level,d_nc,d_ncolumns,d_top_q,d_dp;
    double d_window,d_SDF,d_min_CV,d_min_correl,d_minDE,d_maxDE,d_fudge,d_impute_scale;
    vector<int> d_min_obs,d_ssize,d_ct1,d_ct2,d_np;
    vector<string> d_labels;
    set<string> d_inclusion;


    public:
    Option(const map<string,string>& opm);
    const string& file() const { return d_file; }
    const vector<string>& labels() const { return d_labels; }
    const set<string>& inclusion() const { return d_inclusion; }
    bool rep() const { return d_rep; }
    bool log2() const { return d_log2; }
    bool shared_q() const { return d_shared_q; }
    const string& normalization() const { return d_normalization; }
    double window() const { return d_window; }
    int dp() const { return d_dp; }
    double SDF() const { return d_SDF; }
    double min_CV() const { return d_min_CV; }
    double min_correl() const { return d_min_correl; }
    const vector<int>& min_obs() const { return d_min_obs; }
    const vector<int>& ssize() const { return d_ssize; }
    int min_f() const { return d_min_f; }
    int max_f() const { return d_max_f; }
    int min_q() const { return d_min_q; }
    double minDE() const { return d_minDE; }
    double maxDE() const { return d_maxDE; }
    double fudge() const { return d_fudge; }
    double impute_scale() const { return d_impute_scale; }
    const string& impute_type() const { return d_impute_type; }
    const vector<int>& ct1() const { return d_ct1; }
    const vector<int>& ct2() const { return d_ct2; }
    int nc() const { return d_nc; }
    const vector<int>& np() const { return d_np; }
    int level() const { return d_level; }
    int ncolumns() const { return d_ncolumns; }
    int top_q() const { return d_top_q; }
    int nt() const { return d_ssize.size(); }
};
