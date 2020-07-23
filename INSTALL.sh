#!/usr/bin/bash

base_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA
db_dir=/storage/data/PROJECT/biouser1/TestPaper/obaDIA/db

dir=`pwd`
db=$dir/db

# replace default directories to the current directories

tar zvxf $dir/src/mapDIA.tar.gz >/dev/null
tar zvxf $dir/src/Trinotate-v3.1.0-pro.tar.gz >/dev/null
mv -f mapDIA src/
mv -f Trinotate-v3.1.0-pro src/

rm -rf $dir/bin $dir/env
cp -rf $dir/src/bin $dir/bin
cp -rf $dir/src/env $dir/env

sed -i "s#$base_dir#$dir#g" $dir/bin/oba.pl 
sed -i "s#$base_dir#$dir#g" $dir/env/.bashrc

ls $dir/bin/*/*.*|while read a;do sed -i "s#$base_dir#$dir#g" $a;done
ls $dir/bin/*/*.*|while read a;do sed -i "s#$db_dir#$db#g" $a;done



chmod +x $dir/bin/oba.pl
chmod +x $dir/bin/*/*.*
chmod +x $dir/src/*
chmod +x $dir/src/*/*.*

while [ -n "$1" ]
do
	case "$1" in
		-n) echo "build without signalp" && sed -i 's/.*signalp.*//g' $dir/bin/Anno/annoP.sh && sed -i 's/.*signalp.*//g' $dir/env/.bashrc;;
	esac
	shift
done

echo "finished !"

