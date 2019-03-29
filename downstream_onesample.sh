#!/bin/bash

# Runs filters downstream of Sidron and sorts results into
# groups. Global parameters can be locally overriden by
# writing the new declaration into a file called config.txt
# in the folder where this script is called. A stamp is
# written into a file called log.txt in the same folder.

if [ ! -n "$1" ]
then
  echo "Usage: `basename $0` base_name"
  exit 1
fi

tolog(){
  # Writes something into a file called log.txt
  # If called without arguments, stamps date and time
  if [[ -n $1 && $1 != "" ]]
  then
    echo $1 >>"log.txt"
  else
    echo '----------' >>"log.txt"
    date +"%D %T" >>"log.txt"
  fi
}

echolog(){
  echo $1
  tolog $1
}

bname="$1"

tcov_field=4
tsid_field=8

PSTR=0.00     # Minimum strand probability
CLIMIT=22     # Coverage where sidron's limits change
# S Limit below CLIMIT for both Het and Hz
lHetA=0        # Limit for Het is lHetA*cov+lHetB
lHetB=5.807    #
lHzA=-0.2583   # Limit for Hz is lHzA*cov+lHzB
lHzB=2.6546    #
# S Limit above CLIMIT for both Het and Hz
HetA=0.7019   # Limit for Het is HetA*cov+HetB
HetB=-9.6348  #
HzA=-0.135   # Limit for Hz is HzA*cov+HzB
HzB=-0.056    #
MINCOV=6      # Minimum coverage
MAXHZ=-2      # Maximum score for Hz
# Different zygosities
## One third
MINCOV3=20
HetA3=0.3   # Limit for Het is HetA*cov+HetB
HetB3=-3.369   #
HzA3=0
HzB3=-4


#__________________________________
# End of config-overriden variables
#__________________________________

if [ -e "config.txt" ]
then
  while read fl
  do
  eval $fl
  done <"config.txt"
fi

tolog
tolog "Running downstream.sh $@ with parameters:"
for v in PSTR CLIMIT lHetA lHetB lHzA lHzB HetA HetB HzA HzB MINCOV
do
  eval "val=\$$v"
  tolog "$v=$val"
done

tposs=$(($tsid_field-1))

tcov_var="\$$tcov_field"
tsid_var="\$$tsid_field"
tstr_var="\$$tstr_field"
ts="\$$tposs"

###################
THETSID2=$HetA*$tcov_var+$HetB
THZSID2=$HzA*$tcov_var+$HzB

THETSID1=$lHetA*$tcov_var+$lHetB
THZSID1=$lHzA*$tcov_var+$lHzB

# Awk expressions
passmincov="($tcov_var>=$MINCOV)"
isthet="(($tcov_var<$CLIMIT && $tsid_var>$THETSID1) || ($tcov_var>=$CLIMIT && $tsid_var>$THETSID2))"    # het in tum
isthz="(($tcov_var<$CLIMIT && $tsid_var<$THZSID1) || ($tcov_var>=$CLIMIT && $tsid_var<$THZSID2))"       # hz in tum
                                                                                                          
tnothet="(($tcov_var<$CLIMIT && $tsid_var<$THETSID1) || ($tcov_var>=$CLIMIT && $tsid_var<$THETSID2))"   # tum is not het
eqmfbase="(substr($ts,1,1)==\$3)"                                                           # most frequent base is same as reference



echo "Het variants..."


awk "$passmincov && $isthet" $bname > "$bname.variants"


###################


echo "Hz variants..."
echo "$passmincov && $tnothet && !($eqmfbase)"
awk "$passmincov && $tnothet && !($eqmfbase)" $bname >> "$bname.variants"


