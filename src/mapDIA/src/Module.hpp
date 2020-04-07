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


class Module
{
    bool d_modulebool;
    string d_module_type,d_module,d_mrf;
    int d_min_size,d_max_size,d_module_freq;
    vector<set<int> > d_adj;


    public:
    Module() : d_modulebool(false),
    d_min_size(0),d_max_size(INT_MAX),d_module_freq(0) { }
    void setModulebool(const string& str0);
    bool modulebool() const { return d_modulebool; }
    void setModule(const string& str0);
    const string& module() const { return d_module; }
    void setModule_type(const string& str0);
    const string& module_type() const { return d_module_type; }
    void setMrf(const string& str0);
    const string& mrf() const { return d_mrf; }
    void setModule_size(const string& str0);
    int max_size() const { return d_max_size; }
    int min_size() const { return d_min_size; }
    void setModule_freq(const string& str0);
    int module_freq() const { return d_module_freq; }
    const vector<set<int> >& adj() const { return d_adj; }
    void module_data(const vector<string>& pidvec);
};
