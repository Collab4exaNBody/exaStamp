#!/bin/bash
# Runs a regression test in a fresh copy of its input folder (see AddSingleTest.cmake).
# usage: run_in_test_copy.sh <run_dir> <input_dir> <potentials_dir> <command...>
#  - <run_dir> is deleted and recreated, then filled with a copy of <input_dir> (symlinks followed)
#  - <run_dir>/data points to <potentials_dir>, so that './data' in ONIKA_DATA_PATH finds the potential files
#  - <command...> runs from <run_dir> : outputs stay there, sources are never written
run_dir=$1 ; input_dir=$2 ; potentials_dir=$3
shift 3
case "$run_dir" in
  */regression_runs/*) ;;
  *) echo "run_in_test_copy.sh: refusing to use run directory '$run_dir'" >&2 ; exit 1 ;;
esac
rm -rf "$run_dir" && mkdir -p "$run_dir" || exit 1
cp -rL "$input_dir"/. "$run_dir"/ || { echo "run_in_test_copy.sh: cannot copy '$input_dir'" >&2 ; exit 1 ; }
[ -e "$run_dir/data" ] || ln -s "$potentials_dir" "$run_dir/data"
cd "$run_dir" || exit 1
exec "$@"
