plot_idle()
{
  ( SEP=" plot" ; \
    echo -n "set datafile separator \";\" ;" ; \
    for f in $* ; \
    do \
      echo -n "$SEP \"$f\" every ::1 using 1:2 with lines" ; \
      SEP=" ," ; \
    done \
  ) \
  | xargs -0 -I{} gnuplot -p -e "{}"
}

plot_idle_png()
{
  ( SEP=" plot" ; \
    echo -n "set terminal png transparent truecolor size 1600,1200 linewidth 4 ; set datafile separator \";\" ;" ; \
    for f in $* ; \
    do \
      echo -n "$SEP \"$f\" every ::1 using 1:2 ls -1 with lines" ; \
      SEP=" ," ; \
    done \
  ) \
  | xargs -0 -I{} gnuplot -p -e "{}"
}


