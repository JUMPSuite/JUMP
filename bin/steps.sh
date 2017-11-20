#!/bin/bash
#
this_scripts_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
dtn=`date +%Y%m%d.%N`
pwd="$(pwd)"
this_work_dir="/home/sespy/work/proteomics/jump-search"

echo "start (steps) $this_scripts_dir"
echo "pwd $pwd"
echo "dtn $dtn"

#$this_scripts_dir/step_1.sh $this_work_dir $dtn
#$this_scripts_dir/step_2.sh $this_work_dir $dtn

echo "done. (steps)"


