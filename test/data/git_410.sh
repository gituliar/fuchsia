#/usr/bin/bash

for x1 in "x" "x+1" "x-1";
do
  i1=git_410
  o1=${i1}_$x1
  cmd="delirium maple_super_reduce -o $o1.mtx -p $x1 $i1.mtx"
  echo $cmd
  eval $cmd
  for x2 in "x" "x+1" "x-1";
  do
    if [ $x1 != $x2 ]
    then
      i2=${o1}
      o2=${i2}_$x2
      cmd="delirium maple_super_reduce -o $o2.mtx -p $x2 $i2.mtx"
      echo $cmd
      $cmd
      for x3 in "x" "x+1" "x-1";
      do
        if [[ $x1 != "x" && $x1 != $x3 && $x2 != $x3 ]]
        then
          i3=${o2}
          o3=${i3}_$x3
          cmd="delirium maple_super_reduce -o $o3.mtx -p $x3 $i3.mtx"
          echo $cmd
          $cmd
        fi
      done
    fi
  done
done
