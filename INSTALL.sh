#!/usr/bin/bash

base_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA
db_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db

dir=`pwd`
db=$dir/db

# replace default directories to the current directories

sed -i "s#$base_dir#$dir#g" $dir/bin/oba.pl 
ls $dir/bin/*/*.*|while read a;do sed -i "s#$base_dir#$dir#g" $a;done
ls $dir/bin/*/*.*|while read a;do sed -i "s#$db_dir#$db#g" $a;done

chmod +x $dir/bin/oba.pl
chmod +x $dir/bin/*/*.*
chmod +x $dir/src/*
chmod +x $dir/src/*/*.*

echo "finished !"

