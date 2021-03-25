#!/bin/bash
##########################################
#Merge two csvs by common column
# Created: 08/18/2017 by Feng Xue
##########################################
#Parameters:
#1: file1
#2: file2
#3: common column in file1, 1 by default
#4: common column in file2, 1 by default
#5: output
#6: output format
#7: columns to exclude
#8: Whether to keep differtial data
##########################################
file_1=$1
file_2=$2
common_column_1=$3
common_column_2=$4
mergedfile=$5
output_format=$6
field_remove=$7
nodiff=$8

[ ${#file_1} -eq 0 ] && echo "1st file is empty" && exit -1
[ ${#file_2} -eq 0 ] && echo "2nd file is empty" && exit -1
[ ! -f $file_1 ] && echo "1st file is not readable" && exit -1
[ ! -f $file_2 ] && echo "2nd file is not readable" && exit -1
[ $file_1 = $file_2 ] && "1st file and 2nd file should be different ones" && exit -1

[ ${#common_column_1} -eq 0 ] && common_column_1=1
[ ${#common_column_2} -eq 0 ] && common_column_2=1
int_re='^[0-9]+$'
if ! [[ $common_column_1 =~ $int_re ]] ; then
   echo "common clumns should be integer" ; exit -1
fi
if ! [[ $common_column_2 =~ $int_re ]] ; then
   echo "common clumns should be integer" ; exit -1
fi

[ ${#mergedfile} -eq 0 ] && echo "Please provide filename for the merged file" && exit -1

[ ${#nodiff} -gt 0 ] && nodiff=1 || nodiff=0



field_remove_1=`echo $field_remove| grep -oE "1\.[0-9]*,|1\.[0-9]*$"`
field_remove_2=`echo $field_remove| grep -oE "2\.[0-9]*,|2\.[0-9]*$"`
[ ${#field_remove_1} -eq 0 ] && field_remove_1="1.$common_column_1" || field_remove_1="$field_remove_1,1.$common_column_1"
[ ${#field_remove_2} -eq 0 ] && field_remove_2="2.$common_column_2" || field_remove_2="$field_remove_2,2.$common_column_2"

file_ext=`basename $mergedfile|rev|cut -d. -f 1|rev`
mergedfile_tmp1=`echo $mergedfile|sed -e "s/\.$file_ext$/_tmp1.$file_ext/g"`
mergedfile_tmp2=`echo $mergedfile|sed -e "s/\.$file_ext$/_tmp2.$file_ext/g"`
header_1=`head -n 1 $file_1`
header_2=`head -n 1 $file_2`

output_format_1=`echo $output_format| grep -oE "1\.[0-9]*,|1\.[0-9]*$"`
output_format_2=`echo $output_format| grep -oE "2\.[0-9]*,|2\.[0-9]*$"`
[ ${#output_format_1} -gt 0 ] && has_1=1 || has_1=0
[ ${#output_format_2} -gt 0 ] && has_2=1 || has_2=0

if [ $has_1 -eq 0 ]
then
  column_count_1=`echo $header_1| grep -o "," | wc -l`
  ((column_count_1=column_count_1+1))
  for ((i=1;i<=$column_count_1;i++))
  do
    tmp=`echo $field_remove_1 | grep -E "1\.$i,|1\.$i$"`
    if [ ${#tmp} -eq 0 ]
    then
      if [ $i -lt $column_count_1 ]
      then
        output_format_1="${output_format_1}1.$i,"
      else
        output_format_1="${output_format_1}1.$i"
      fi
    fi
  done
  output_format_1=`echo $output_format_1|sed -e "s/,$//g"`
  output_format="$output_format_1,$output_format"
fi
if [ $has_2 -eq 0 ]
then
  column_count_2=`echo $header_2 | grep -o "," | wc -l`
  ((column_count_2=column_count_2+1))
  for ((i=1;i<=$column_count_2;i++))
  do
    tmp=`echo $field_remove_2 | grep -E "2\.$i,|2\.$i$"`
    if [ ${#tmp} -eq 0 ]
    then
      if [ $i -lt $column_count_2 ]
      then
        output_format_2="${output_format_2}2.$i,"
      else
        output_format_2="${output_format_2}2.$i"
      fi
    fi
  done
  output_format_2=`echo $output_format_2|sed -e "s/,$//g"`
  output_format="$output_format,$output_format_2"
fi
output_format=`echo $output_format|sed -e "s/,,/,/g" -e "s/ //g" -e "s/,$//g"`
join --header -t',' -1 $common_column_1 -2 $common_column_2 -a1 -o 0,$output_format <( echo $header_1 && tail -n +2 $file_1 |sort -t, -k$common_column_1) <(echo $header_2 && tail -n +2 $file_2 | sort -t, -k$common_column_2) > $mergedfile_tmp1

output_format=`echo $output_format|sed -e "s/1\./3./g" -e "s/2\./1./g"`
output_format=`echo $output_format|sed -e "s/3\./2./g"`
join --header -t',' -1 $common_column_2 -2 $common_column_1 -a1 -o 0,$output_format <( echo $header_2 && tail -n +2 $file_2 |sort -t, -k$common_column_2) <(echo $header_1 && tail -n +2 $file_1 | sort -t, -k$common_column_1) > $mergedfile_tmp2

grep -axFf $mergedfile_tmp1 $mergedfile_tmp2 >$mergedfile

if [ $nodiff -eq 0 ]
then
  grep -avxFf $mergedfile_tmp1 $mergedfile_tmp2 >>$mergedfile
  grep -avxFf $mergedfile_tmp2 $mergedfile_tmp1 >>$mergedfile
fi
rm -f $mergedfile_tmp1 $mergedfile_tmp2
