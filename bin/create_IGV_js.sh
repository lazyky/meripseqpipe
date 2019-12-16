#! /bin/bash
fasta=$1
gtf=$2
merged_peak_file=$3
designfile=$4
echo "Start to generate IGV.js"

## setting tmp files' name
bigwig_tracks_file=tmp.bigwig.tracks
peaks_tracks_file=tmp.peaks.tracks

## combined tracks of bigwig
sampleinfo_list=$(awk 'BEGIN{FS=","}NR>1{print $1","$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
for sample_group_id in ${sampleinfo_list}
    do
    {  
        sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
        group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
        bigwig_input_file=$(ls ${sample_id}.input_*.igv.bigwig)
        bigwig_ip_file=$(ls ${sample_id}.ip_*.igv.bigwig)
        cat >> ${bigwig_tracks_file} << EOF
        {
            url: 'http://localhost:8080/${bigwig_input_file}',
            name: '${sample_id}.input',
            color: 'rgb(200,0,0)',
            autoscaleGroup: 'group_${group_id}.${sample_id}'
        },
        {
            url: 'http://localhost:8080/${bigwig_ip_file}',
            name: '${sample_id}.ip',
            color: 'rgb(200,0,0)',
            autoscaleGroup: 'group_${group_id}.${sample_id}'
        },
EOF
    }
done

## combined tracks of merged group peaks
groups_peak_file=$(ls *_merged_group_*igv.bed)
for peak_file in ${groups_peak_file}
    do
    {  
        cat >> ${peaks_tracks_file} << EOF
        {
            type: "annotation",
            format: "bed",
            url: 'http://localhost:8080/${peak_file}',
            name: "${peak_file}"
        },
EOF
    }
done

## combined tracks and allpeaks track
cat ${bigwig_tracks_file} ${peaks_tracks_file} > tmp.tracks
cat >> tmp.tracks << EOF
        {
            type: "annotation",
            format: "bed",
            url: 'http://localhost:8080/${merged_peak_file}',
            name: "${merged_peak_file}"
        }
EOF
tracks_js=$(cat tmp.tracks)

## combined all info 
cat>igv.js<<EOF
var igvDiv = document.getElementById("igvDiv");
var options =
{
    reference: {
        id: "$fasta",
        fastaURL: "http://localhost:8080/hg38_chr22.fa",
        indexURL: "http://localhost:8080/hg38_chr22.fa.fai",
        wholeGenomeView: false
    },
    locus: 'chr22',
    tracks: [
        {
            type: "annotation",
            format: "gtf",
            url: "http://localhost:8080/$gtf",
            displayMode: "SQUISHED",
            name: "$gtf",
            visibilityWindow: 10000000
        },
$tracks_js
]
}; 
var browser = igv.createBrowser(igvDiv, options);
EOF
rm -rf tmp*
echo "Generate IGV.js was success"