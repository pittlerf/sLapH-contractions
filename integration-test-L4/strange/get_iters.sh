echo "# $(date)" > iters.dat
echo "conf rnd iters" >> iters.dat
for cnf in cnfg*; do
  conf=$( echo $cnf | sed "s/cnfg//g" )
  echo ${cnf}
  for rnd in 00 01 02 03 04 05; do
    output=${cnf}/outputs/*${rnd}.out
    if [ -f ${output} ]; then
      iters=$(grep Done ${output} | grep iter | awk '{print $4}')
      for i in $iters; do
        echo $conf $rnd ${i} >> iters.dat
      done
    fi
  done
done
