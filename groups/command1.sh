IFS=$'\n'; for i in `ls -1 *.txt`; do Rscript do.one.R "$i"; done
