#!/bin/bash
#
#this_work_dir=$1
this_scripts_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
dtn=`date +%Y%m%d.%N`
pwd="$(pwd)"
this_work_dir="/home/sespy/work/proteomics/jump-search"

echo "start (steps) $this_scripts_dir"
echo "pwd $pwd"
echo "this_scripts_dir $this_scripts_dir"
echo "this_work_dir $this_work_dir"
echo "dtn $dtn"

$this_work_dir/bin/step_1.sh $this_work_dir $dtn
$this_work_dir/bin/step_2.sh $this_work_dir $dtn

echo "done. (steps)"


