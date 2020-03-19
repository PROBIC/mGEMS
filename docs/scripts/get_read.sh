ftppath="ftp://ftp.sra.ebi.ac.uk/vol1/fastq"

strlen=$((${#1}))
dir1="/"${1:0:6}

if [ ${strlen} -lt 10 ]
then
dir2="/"
elif [ ${strlen} -lt 11 ]
then
dir2="/00"${1: -1}"/"
elif [ ${strlen} -lt 12 ]
then
dir2="/0"${1: -2}"/"
elif [ ${strlen} -lt 13 ]
then
dir2="/"${1: -3}"/"
else
echo "check accession number"
fi

dlpath=$ftppath$dir1$dir2$1"/"
echo "downloading reads from path: "$dlpath

wget $dlpath"*"
