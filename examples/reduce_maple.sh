#/usr/bin/bash

for x1 in "x" "x+1" "x-1";
do
  i1=$1
  o1=${i1}_maple_$x1
  cmd="fuchsia reduce --maple -o $o1.mtx -p $x1 $i1.mtx"
  echo $cmd
  time $cmd
  for x2 in "x" "x+1" "x-1";
  do
    if [ $x2 != $x1 ]
    then
      i2=${o1}
      o2=${i2}_$x2
      cmd="fuchsia reduce --maple -o $o2.mtx -p $x2 $i2.mtx"
      echo $cmd
      time $cmd
      for x3 in "x" "x+1" "x-1";
      do
        if [[ $x3 != $x1 && $x3 != $x2 ]]
        then
          i3=${o2}
          o3=${i3}_$x3
          cmd="fuchsia reduce --maple -o $o3.mtx -p $x3 $i3.mtx"
          echo $cmd
          time $cmd
        fi
      done
    fi
  done
done
