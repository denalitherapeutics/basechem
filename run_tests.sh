#!/bin/sh
help()
{
    # Display Help
    echo "Run the test suite and generate the coverage report."
    echo
    echo "usage: ./run_tests.sh [-h] [--exclude-external-pkgs]"
    echo
    echo "options:"
    echo " -h, --help                      show this help message and exit"
    echo " -noi, --exclude-inductive-bio   exclude tests that rely on InductiveBio"
    echo " -nos, --exclude-slurm           exclude tests that rely on Slurm"
    echo " -nod, --exclude-dtx             exclude tests that rely on Dotmatics"
    echo " -e, --exclude-external-pkgs     run only the tests that don't require external packages (dtxwrapper, InductiveBio)"
}
# Function to handle options and arguments
handle_options() {
  while [ $# -gt 0 ]; do
    case $1 in
        -h | --help)
            help
            exit 0
            ;;
        -e | --exclude-external-pkgs )
            excludeArgs+=("--exclude-tag=external")
            ;;
        -noi | --exclude-inductive-bio )
            excludeArgs+=("--exclude-tag=inductive")
            ;;
        -nos | --exclude-slurm )
            excludeArgs+=("--exclude-tag=slurm")
            ;;
        -nod | --exclude-dtx )
            excludeArgs+=("--exclude-tag=dtx")
            ;;
        *)
            echo
            echo "Invalid option: $1" >&2
            echo
            help
            exit 1
            ;;
        esac
    shift
  done
}

###### Begin script execution ######
excludeArgs=()
handle_options "$@"

set -ex

black --exclude='roles|node_modules|\.git|.*cache.*|\.tox|venv|common/gypsum_dl|common/apbs|common/ESP_DNN' --check . 
isort . --profile black --skip-glob "**/common/gypsum_dl" --skip-glob "**/common/apbs" --skip-glob "**/common/ESP_DNN" --check-only

coverage erase
python manage.py collectstatic --noinput # needed for mni_common tests
python manage.py makemigrations --dry-run --check

coverage run manage.py test --keepdb --debug-mode "${excludeArgs[@]}" --noinput
coverage report -m --skip-covered --fail-under 85 --omit="basechem/common/gypsum_dl/*,basechem/common/ESP_DNN/*,basechem/main/admin*.py,basechem/common/mmpdb_utils.py" > "basechem/common/tests/coverage_report.txt"
