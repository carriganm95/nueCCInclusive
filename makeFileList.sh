#!/bin/bash

htgettoken -a htvaultprod.fnal.gov -i icarus

definition="micarrig_shower_dEdx_cal_reco1_to_caf_shower_dEdx_cal"
#definition="Icaruspro_2025_Run2_production_Run2reprocess_v09_89_01_02p02_offbeamnumimajority_flatcaf_prescaled" #cosmic data
localPath="/pnfs/sbn/data/sbn_fd/poms_production/2025A_icarus_NuMI_MC/FHC_NuMI/mc/reconstructed/icaruscode_v09_89_01_02p02/flatcaf/*/*/*.flat.caf*.root" #mc data
outputList="mcFiles.list"
treeName="recTree"
checkFile=false
debug=false

declare -A existing_files=()

getToken() {
    htgettoken -a htvaultprod.fnal.gov -i icarus
}

# check for existing files in output list
load_existing_files() {
    if [[ -f "$outputList" ]]; then
        while IFS= read -r line; do
            [[ -n "$line" ]] && existing_files["$line"]=1
        done < "$outputList"
    else
        touch "$outputList"
    fi
}

# function to check if a ROOT file can be opened and contains a non-empty tree
checkRootFile() {
    local filePath="$1"
    local treeName="$2"
    
    if ! root -l -b -q -e "TFile* tf = TFile::Open(\"$filePath\"); if (!tf || tf->IsZombie()) gSystem->Exit(1); TTree* tree = nullptr; tf->GetObject(\"$treeName\", tree); if (!tree) gSystem->Exit(2); if (tree->GetEntries() == 0) gSystem->Exit(3);" &>/dev/null; then
        status=$?
        case $status in
            1) echo "Warning: failed to open $filePath with ROOT" >&2 ;;
            2) echo "Warning: tree '$treeName' not found in $filePath" >&2 ;;
            3) echo "Warning: tree '$treeName' empty in $filePath" >&2 ;;
            *) echo "Warning: ROOT check failed for $filePath (code $status)" >&2 ;;
        esac
        return 1
    fi
    return 0
}

# function to process local files
process_files_local() {

    echo "Processing local files from $localPath"

    shopt -s nullglob

    local file_pattern="${localPath##*/}"
    local prefix_before_glob="${localPath%%\**}"
    local remainder="${localPath#${prefix_before_glob}}"
    local first_segment="${remainder%%/*}"
    local first_layer_glob

    if [[ "${remainder}" == "${localPath}" ]]; then
        # No wildcard found; fall back to the directory containing the files
        first_layer_glob="${localPath%/*}"
    else
        if [[ -z "${first_segment}" || "${first_segment}" == "${remainder}" ]]; then
            first_layer_glob="${localPath%/*}"
        else
            first_layer_glob="${prefix_before_glob}${first_segment}"
        fi
    fi

    local counter=0
    local total_dirs=0
    local total_candidates=0
    local stop_processing=false

    for dir in ${first_layer_glob}; do
        if [[ ! -d "$dir" ]]; then
            continue
        fi

        total_dirs=$((total_dirs + 1))
        echo "Scanning directory: $dir"

        local dir_candidate_count=0

        while IFS= read -r -d '' f; do
            dir_candidate_count=$((dir_candidate_count + 1))

            if $debug && [ $counter -ge 100 ]; then
                echo "Debug mode active: stopping after $counter files"
                stop_processing=true
                break
            fi

            if [[ ! -f "$f" ]]; then
                echo "Warning: $f not found" >&2
                continue
            fi

            # if [ $counter % 100 -eq 0 ]; then
            #     getToken
            # fi

            local remote_path="$f"
            if [[ "$f" == /pnfs/* ]]; then
                remote_path="root://fndcadoor.fnal.gov://${f#/pnfs/}"
            fi

            if [[ -n "${existing_files["$remote_path"]}" ]]; then
                if $debug; then
                    echo "Skipping existing file: $remote_path"
                fi
                continue
            fi

            if [ $checkFile = true ]; then
                if ! checkRootFile "$remote_path" "$treeName"; then
                    continue
                fi
            fi

            if $debug; then
                echo "Adding new file: $remote_path"
            fi

            echo "$remote_path" >> "$outputList"
            existing_files["$remote_path"]=1
            counter=$((counter + 1))
        done < <(find "$dir" -type f -name "$file_pattern" -print0)

        if (( dir_candidate_count == 0 )); then
            echo "  No files matching ${file_pattern}"
        fi

        total_candidates=$((total_candidates + dir_candidate_count))

        if $stop_processing; then
            break
        fi
    done

    shopt -u nullglob

    echo "Scanned $total_dirs directories (${total_candidates} candidate files)"
    echo "Total new files added: $counter"
}

# function to process SAM files
process_files_sam() {

    files=$(samweb list-definition-files "$definition")

    local counter=0

    if [ counter % 100 -eq 0 ]; then
        getToken
    fi

    for f in $files; do
        if [ $counter -ge 10 ]; then
            break
        fi

        loc=$(samweb locate-file "$f" | awk '{print $1; exit}')
        if [[ -z "$loc" ]]; then
            echo "Warning: no location found for $f" >&2
            continue
        fi

        root_uri=${loc#dcache:/pnfs}
        root_uri="root://fndcadoor.fnal.gov:${root_uri}"
        full_uri="${root_uri}/${f}"

        if [[ -n "${existing_files["$full_uri"]}" ]]; then
            if $debug; then
                echo "Skipping existing file: $full_uri"
            fi
            continue
        fi

        if ! checkRootFile "$full_uri" "$treeName"; then
            continue
        fi

        if $debug; then
            echo "Adding new file: $full_uri"
        fi

        echo "$full_uri" >> "$outputList"
        existing_files["$full_uri"]=1
        counter=$((counter + 1))
    done

    echo "Total new files added: $counter"
}

load_existing_files
process_files_local
#process_files_sam