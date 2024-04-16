IFS=$'\n'; for i in `ls -1`; do echo $i; cat "$i"|sed "1d"; done|paste - -|sort -t" " -g -k2 -r|less
