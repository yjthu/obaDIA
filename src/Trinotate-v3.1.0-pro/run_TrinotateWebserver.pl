#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

my $usage = "usage: $0 [port_no=8080]\n\n";

my $port_no = $ARGV[0] or die $usage;


## if the port comes up as already being used, check what's using it like so:
## sudo fuser -v 3000/tcp
##  and kill it or choose a different port to use.
##
## other useful things to know:
##  check for open ports:
##        sudo netstat -tulpn



 main: {

     my $lighttpd_prog = `which lighttpd`;
     chomp $lighttpd_prog;
     unless ($lighttpd_prog =~ /\w/) {
         die "Error, cannot locate 'lighttpd' program. Be sure to have it installed and accessible from your PATH env var";
     }
     
     my $perl_path = `which perl`;
     chomp $perl_path;
     unless ($perl_path =~ /^\//) {
         die "Error, can't determine where to find perl....  'which perl' returns: $perl_path ";
     }
     
     my $document_root = "$FindBin::RealBin/TrinotateWeb";
     my $conf_file_template = "$FindBin::RealBin/TrinotateWeb.conf/lighttpd.conf.template";
     
     my $conf_file = "$FindBin::RealBin/TrinotateWeb.conf/lighttpd.conf.port$port_no";
     unless (-s $conf_file) {
         my $template = `cat $conf_file_template`;
         $template =~ s/__DOCUMENT_ROOT__/$document_root/ or die "Error, could not replace __DOCUMENT_ROOT__ in $conf_file_template";
         $template =~ s/__PORT_NO__/$port_no/ or die "Error, could not replace __PORT_NO__ in $conf_file_template";

         $template =~ s/__PERL_PATH__/$perl_path/g or die "Error, could not replace __PERL_PATH__ in $conf_file_template with $perl_path";
         
         open (my $ofh, ">$conf_file") or die "Error, cannot write to $conf_file";
         print $ofh $template;
         close $ofh;
     }

     my $cmd = "$lighttpd_prog -D -f $conf_file";
     print STDERR "$cmd\n";
     my $ret = system($cmd);

     exit($ret);
}

