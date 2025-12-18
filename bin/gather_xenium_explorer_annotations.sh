#!/usr/bin/env bash

#################
### CONSTANTS ###
#################

SELECTION_NAME_STR="Selection name"

##################
### INPUT ARGS ###
##################

input_files=""
output_file="manual_annotations.tsv"

while [[ $# -gt 0 ]]; do
    flag=$1

    case "${flag}" in
        -i) input_files=$input_files" "$2; shift;;
        -o) output_file=$2; shift;;
        *) echo "Unknown option $1 ${reset}" && exit 1
    esac

    shift
done

if [ -z "${input_files}" ]; then
    echo "No input files provided. Use -i to specify input files."
    exit 1
fi

###############################
### CREATE ANNOTATIONS FILE ###
###############################

rm -f "${output_file}"
echo -e "Cell_ID\tTissue_annotation" > "${output_file}"

echo $input_files

echo $input_files | while read input_file; do
    annotation_class=$(grep "$SELECTION_NAME_STR" "${input_file}" | cut -f2 -d ':' | tr -d ' ')

    egrep -v '^#|Cell ID' "${input_file}" | \
        cut -f1 -d',' | \
        sed "s/\$/\t${annotation_class}/g" >> "${output_file}"
done
