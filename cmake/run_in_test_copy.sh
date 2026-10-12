#!/bin/bash
# Runs a regression test in the run directory of its input folder (see AddSingleTest.cmake).
# usage: run_in_test_copy.sh <run_dir> <input_dir> <potentials_dir> <command...>
#  - <run_dir> is shared by all tests of <input_dir> and kept between runs : some tests read files written by
#    another test of the folder (generate/read pairs), check_values compares against the reference file written
#    by the first run
#  - the content of <input_dir> is copied into <run_dir> (symlinks followed, files refreshed when the source is newer)
#  - <run_dir>/data points to <potentials_dir>, so that './data' in ONIKA_DATA_PATH finds the potential files
#  - <command...> runs from <run_dir> : sources are never written
run_dir=$1 ; input_dir=$2 ; potentials_dir=$3
shift 3
case "$run_dir" in
  */regression_runs/*) ;;
  *) echo "run_in_test_copy.sh: refusing to use run directory '$run_dir'" >&2 ; exit 1 ;;
esac
mkdir -p "$run_dir" || exit 1
copy_inputs()
{
  cp -rLu "$input_dir"/. "$run_dir"/ || return 1
  [ -e "$run_dir/data" ] || ln -s "$potentials_dir" "$run_dir/data"
}
if command -v flock > /dev/null ; then
  # tests of the same folder may start together (ctest -j)
  ( flock 9 && copy_inputs ) 9> "$run_dir/.copy.lock"
else
  copy_inputs
fi || { echo "run_in_test_copy.sh: cannot copy '$input_dir' to '$run_dir'" >&2 ; exit 1 ; }
cd "$run_dir" || exit 1
exec "$@"
