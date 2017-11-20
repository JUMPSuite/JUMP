#!/bin/bash
#
this_scripts_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
dtn=`date +%Y%m%d.%N`
pwd="$(pwd)"

#

echo "start (steps)" $this_scripts_dir
echo "pwd" $pwd

$this_scripts_dir/step_1.sh $dtn
$this_scripts_dir/step_2.sh $dtn

echo "done. (steps)"


