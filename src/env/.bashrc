# .bashrc

# User specific aliases and functions

alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
alias cpan=/usr/local/bin/cpan
# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

# obaDIA enviroment
# set 'signalp_home' correctly and 
# copy the following content to ~/.bashrc file

oba_home=/storage/data/PROJECT/biouser1/TestPaper/obaDIA
export PATH=$oba_home/src:$oba_home/src/mapDIA:$oba_home/src/Trinotate-v3.1.0-pro:$PATH

signalp_home=/storage/data/PUBLIC/softwares/SignalP/signalp-4.1
export PATH=:$signalp_home:$PATH
export PERL5LIB=$signalp_home/lib/:$PERL5LIB



