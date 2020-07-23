#!/usr/bin/bash

fa=$1
out=$2

awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0" ":$0 }'  $1 |awk '{print $1"\t"length($NF)}'|sed 's/>//' >$2
