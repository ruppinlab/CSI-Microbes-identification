#######################################################################################################################
#Autor: Sandro Valenzuela (sandrolvalenzuelad@gmail.com)
#pasteTaxID have 2 usage
#Usage 1: bash parseTaxID.bash --workdir [fastas_path] if you have a lot fastas in the workdir
#Usage 2: bash parseTaxID.bash --multifasta [multifasta_file] if you have a huge multifasta file (.fna, .fn works too)
#for debug add --debug
#######################################################################################################################

if [[ "$@" =~ "--debug" ]]; then
	set -ex
else
	set -e
fi

function makePythonWork {
	echo '#!/usr/bin/python

import glob,os
import sys

if len(sys.argv) >= 3:
	os.chdir(sys.argv[1])
	fastaapend=open(sys.argv[2],"wb")

	for file in glob.glob("*.fasta"):
		fasta=open(file, "r")
		oneline = fasta.readline()
		seq =(file,oneline)
		s=" "
		fastaapend.write(s.join(seq))
		fasta.close()

	fastaapend.close()

else:
        print "Workpath and file_out are needed";' > appendheaders.py
}

function makeAwkWork {
	echo 'BEGIN{FS="|"}
{
if($1~">"){
	split($1,array,">")
	$1=array[2]
	for (i=1;i<100;i++){
	band=0;
		if($i != ""){
				if($i==ID){
					ID=$(i+1);
					print ID
					exit 0
				}
		}

	}
}
}' > parsefasta.awk
}

function makeMergeWork {
echo '#!/usr/bin/python

import glob,os
import sys

if len(sys.argv) >= 3:
	os.chdir(sys.argv[1])
	fastappend=open(sys.argv[2], "w")
	for file in glob.glob("*.fasta"):
		fasta=open(file, "r")
		lines = fasta.readlines()
		seq="".join(lines)
		fastappend.write(seq)
		fasta.close()
	fastappend.close()
else:
	print "Workpath and file_out are needed";' > merge.py
}

function fetchFunction () {
echo '
if [[ "$@" =~ "--debug" ]]; then
	set -ex
else
	set -e
fi
apkikey=$4
function notFoundMessage () {
	if [ "$1" == "" ];then
		echo "###### Warning: no ti found for $2 ######"
	fi
}
headers=$1
switchfile=$2
total=$(wc -l $headers |awk '\''{print $1}'\'')
declare pids
i=1
connectionRetries=2
cat $headers |while read line
do
	echo "* process $3: fetching taxid $i of $total"
	#first, we get the critical data through awk and the ID that we find
	fasta=$(echo $line |awk '\''{print $1}'\'')
	fastaheader=$(echo $line |awk '\''{print $2}'\'')
	ti=$(echo "$fastaheader" |awk -v ID="ti" -f parsefasta.awk)
	opcion="ti"
	if [ "$ti" == "" ];then
		acc=$(echo "$fastaheader" |awk -v ID="acc" -f parsefasta.awk)
		opcion="acc"
			if [ "$acc" == "" ];then
					gi=$(echo "$fastaheader" |awk -v ID="gi" -f parsefasta.awk)
					opcion="gi"
					if [ "$gi" == "" ];then
							gb=$(echo "$fastaheader" |awk -v ID="gb" -f parsefasta.awk)
							opcion="gb"
							if [ "$gb" == "" ];then
									emb=$(echo "$fastaheader" |awk -v ID="emb" -f parsefasta.awk)
									opcion="emb"
									if [ "$emb" == "" ];then
										ref=$(echo "$fastaheader" |awk -v ID="ref" -f parsefasta.awk)
										opcion="ref"
										if [ "$ref" == "" ];then
											opcion=""
										fi
									fi
							fi
					fi
			fi
	fi
	case $opcion in
		"ti")
			echo "* Tax Id already found in $fasta, continue"
			echo "$fasta $ti" >> $switchfile
		;;
		"acc")
			retry=$connectionRetries
			while [ "$ti" == "" ] | [ $retry -ge 1 ]
			do
				ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&${apkikey}&id=$acc" | grep "<Id>"| tail -n1 | awk '\''{print $1}'\'' | cut -d '\''>'\'' -f 2 | cut -d '\''<'\'' -f 1)
				if [ "$ti" == "" ];then
					ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?dbfrom=taxonomy&${apkikey}&id=$acc&rettype=fasta&retmode=xml" | head -n10 | grep "TSeq_taxid" | cut -d '\''>'\'' -f 2 | cut -d '\''<'\'' -f 1 )
				fi
				retry=$((retry-1))
			done
			notFoundMessage $ti $acc
			echo "$fasta $ti" >> $switchfile
		;;
		"gi")
			retry=$connectionRetries
			while [ "$ti" == "" ] | [ $retry -ge 1 ]
			do
				ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&${apkikey}&id=$gi" |grep "<Id>"|tail -n1 |awk '\''{print $1}'\'' |cut -d '\''>'\'' -f 2 |cut -d '\''<'\'' -f 1)

				if [ "$ti" == "$gi" ];then
					ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&${apkikey}&id=$gi" |head -n20 |grep "id" |awk '\''{print $2}'\'' |head -n1)
				fi
				retry=$((retry-1))
			done
			notFoundMessage $ti $gi
			echo "$fasta $ti" >> $switchfile
		;;
		"gb")
			retry=$connectionRetries
			while [ "$ti" == "" ] | [ $retry -ge 1 ]
			do
				gi=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&${apkikey}&id=$gb&rettype=fasta" |awk -v ID="gi" -f parsefasta.awk)
				if [ "$gi" == "" ];then
					ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&${apkikey}&id=$gb" |grep "<Id>"|tail -n1 |awk '\''{print $1}'\'' |cut -d '\''>'\'' -f 2 |cut -d '\''<'\'' -f 1)
				else
					ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&${apkikey}&id=$gi" |grep "<Id>"|tail -n1 |awk '\''{print $1}'\'' |cut -d '\''>'\'' -f 2 |cut -d '\''<'\'' -f 1)
				fi
				retry=$((retry-1))
			done
			notFoundMessage $ti $gb
			echo "$fasta $ti" >> $switchfile
		;;
		"emb")
			retry=$connectionRetries
			while [ "$ti" == "" ] | [ $retry -ge 1 ]
			do
				gi=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&${apkikey}&id=$emb&rettype=fasta" |awk -v ID="gi" -f parsefasta.awk)
				ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=taxonomy&${apkikey}&id=$gi" |grep "<Id>"|tail -n1 |awk '\''{print $1}'\'' |cut -d '\''>'\'' -f 2 |cut -d '\''<'\'' -f 1)
				retry=$((retry-1))
			done
			notFoundMessage $ti $emb
			echo "$fasta $ti" >> $switchfile
		;;
		"ref")
			retry=$connectionRetries
			while [ "$ti" == "" ] | [ $retry -ge 1 ]
			do
				ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&${apkikey}&id=$ref&rettype=fasta&retmode=xml" |grep "TSeq_taxid" |cut -d '\''>'\'' -f 2 |cut -d '\''<'\'' -f 1 )
				retry=$((retry-1))
			done
			notFoundMessage $ti $ref
			echo "$fasta $ti" >> $switchfile
		;;
		*)
			#in case the id is not found
			#trying first string as Accession number

			retry=$connectionRetries
			ac=$(echo "$fastaheader" |awk '\''{gsub(">","");print $1}'\'')
			ti=""
			while [ "$ti" == "" ] && [ $retry -ge 1 ]
			do
				ti=$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sequences&${apkikey}&id=$ac&rettype=fasta&retmode=xml" |head -n10 |grep "TSeq_taxid" | cut -d">" -f 2 | cut -d"<" -f 1)
				retry=$((retry-1))
			done

			if [  "$ti" != "" ];then
				echo "$fasta $ti" >> $switchfile
			else
				echo "No id to fetch is available for $ac, continue"
			fi
		;;
	esac


	i=$((i+1))
done' > fetch.bash

}

statusband=0
workpathband=0
multifband=0
multiway=0
pbin=0
PYTHONBIN=/usr/bin/python
parallelJ=1
parallelband=0
apikeyband=0
apkikey=""
unameOut="$(uname -s)"

for i in "$@"
do
	case $i in
	"--workdir")
		workpathband=1
	;;
	"--multifasta")
		multifband=1
	;;
	"--pythonBin")
		phome=1
	;;
	"--parallelJobs")
		parallelband=1
	;;
	"--apikey")
		apikeyband=1
	;;
	"--help")
		echo "UPDATE FROM NCBI: SINCE DECEMBER 2018, NO MORE THAN THREE PARALLEL JOBS CAN FETCH THE DBs."
		echo -e "FOR MORE THAN THREE JOBS PLEASE CREATE AN NCBI ACOUNT AND THEN AN API KEY (USE THE --apikey)\n how to create it?: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/\n\n"
		echo "Usage 1: bash parseTaxID.bash --workdir [fastas_path] if you have a lot fastas in the workdir"
		echo "Usage 2: bash parseTaxID.bash --multifasta [multifasta_file] if you have a huge multifasta file (.fna, .fn works too)"
		echo "Usage 3: bash parseTaxID.bash --multifasta [multifasta_file] --pythonBin to provide a python v2.7"
		echo "Usage 4: bash parseTaxID.bash --multifasta [multifasta_file] --parallelJobs 10 to fetch 10 tax IDs at the same time (default 1, max [all cores or 3 without api key])"
		echo "Usage 5: bash parseTaxID.bash --workdir [fastas_path] --apkikey 222220a2f4875347e41dfc568 (just an example of api key)"
		exit

	;;
	"-h")
		echo "UPDATE FROM NCBI: SINCE DECEMBER 2018, NO MORE THAN THREE PARALLEL JOBS CAN FETCH THE DBs."
		echo -e "FOR MORE THAN THREE JOBS PLEASE CREATE AN NCBI ACOUNT AND THEN AN API KEY (USE THE --apikey)\n how to create it?: https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/\n\n"
		echo "Usage 1: bash parseTaxID.bash --workdir [fastas_path] if you have a lot fastas in the workdir"
		echo "Usage 2: bash parseTaxID.bash --multifasta [multifasta_file] if you have a huge multifasta file (.fna, .fn works too)"
		echo "Usage 3: bash parseTaxID.bash --multifasta [multifasta_file] --pythonBin to provide a python v2.7"
		echo "Usage 4: bash parseTaxID.bash --multifasta [multifasta_file] --parallelJobs 10 to fetch 10 tax IDs at the same time (default 1, max [all cores or 3 without api key])"
		echo "Usage 5: bash parseTaxID.bash --workdir [fastas_path] --apkikey 222220a2f4875347e41dfc568 (just an example of api key)"
		exit

	;;
	*)

		if [ $((workpathband)) -eq 1 ];then
			statusband=$((statusband+1))
			workpathband=0
			WORKDIR=$i
			EXECUTEWORKDIR=$(pwd)
		fi

		if [ $((multifband)) -eq 1 ];then
			statusband=$((statusband+1))
			multifband=0
			multiway=1
			#actual=$(pwd)
			multifname=$(echo "$i" |rev |cut -d '/' -f 1 |rev)
			#multffolder=$(echo "$i" |rev |cut -d '/' -f 2- |rev)
			#if [ "$multifname" == "$multffolder" ];then
			#
			#else
			#	cd $multffolder
			#	multffolder=$(pwd)
			#	cd $actual
			#	multif="$multffolder/$multifname"
			#fi
			multif=$(realpath $i)
			if [ ! -f "$multif" ];then
				echo "* Error: $multif doesn't exist"
				exit
			fi
		fi

		if [ $((pbin)) -eq 1 ];then
			phome=0
			PYTHONBIN=$i
		fi

		if [ $((parallelband)) -eq 1 ];then
			parallelband=0
			parallelJ=$i
			if [ $((parallelJ)) -le 0 ];then
				parallelJ=1
			fi
			if [ $((parallelJ)) -ge $(nproc) ];then
				echo "* Warning: parallelJobs values will set to max CPUs detected ($(nproc))"
				parallelJ=$(nproc)
			fi
		fi

		if [ $((apikeyband)) -eq 1 ];then
			apikeyband=0
			apkikey=$(echo "api_key="$i"&")
		fi

	esac
done

if [ $((statusband)) -eq 1 ]; then

######################		SPLIT FASTAS	##########################
	case $multiway in
	"0")
		echo "* working with workdir parameter"
	;;
	"1")
		echo "* Splitting multifasta, (if the file is a huge file (~300.000 or more sequences), you should go for a coffee while the script works"
		if [ -f $multif ];then
			rm -fr $multifname""_TMP_FOLDER_DONT_TOUCH
			mkdir $multifname""_TMP_FOLDER_DONT_TOUCH
			cd $multifname""_TMP_FOLDER_DONT_TOUCH
			awk '/^>/{close(s);s=++d".fasta"} {print > s}' $multif
			echo "* Splitting complete, DON'T TOUCH $multifname_TMP_FOLDER WHILE SCRIPT IS RUNNING"
			WORKDIR=$(pwd)
			cd ..

		else
			echo "exist($multifasta) = FALSE"
			exit
		fi

	;;
	*)
		echo "unknow error at multifastaflag"
		exit
	;;
	esac


######################		MAKE HEADERS	##########################


	fileout="headers.txt"
	echo "* Making headers from fastas"
	cd $WORKDIR
	makePythonWork
	case $multiway in
		"0")
			#workpath
			$PYTHONBIN appendheaders.py "." $fileout	#just take the first line of each fasta (>foo|1234|lorem ipsum)
		;;
		"1")
			#multif
			$PYTHONBIN appendheaders.py $WORKDIR $fileout	#just take the first line of each fasta (>foo|1234|lorem ipsum)
		;;
	esac

######################################################################
######################		FETCH ID		##########################

	switchfile="newheader.txt"
	touch $switchfile
	total=$(wc -l $fileout |awk '{print $1}')
	makeAwkWork
	fetchFunction

	if [ $((parallelJ)) -ge  2 ]; then
		if [ "$apikey" == "" ] && [ $((parallelJ)) -eq  2 ] ;then
			parallelJ=2
		fi
		total=$(echo $total |awk -v parallelJ=$parallelJ '{print int($1/parallelJ)+1}' )
		split -l $total $fileout
		declare gpids
		i=0
		for Xchunks in $(ls -1 x[a-z][a-z])
		do
			if [[ "$@" =~ "--debug" ]]; then
				bash fetch.bash $Xchunks $switchfile $i $apkikey "--debug" & gpids[${i}]=$(echo $!)
			else
				bash fetch.bash $Xchunks $switchfile $i $apkikey & gpids[${i}]=$(echo $!)
			fi

			i=$((i+1))
		done

		for id in ${gpids[@]}
		do
			case "${unameOut}" in
			    Linux*)
					   wait $id
				;;
			    Darwin*)
					while kill -0 $id >/dev/null 2>&1
					do
					    sleep 10
					done
				;;
			    *)
					echo "Not compatible OS"
			esac

		done
		rm x[a-z][a-z]

	else
		if [[ "$@" =~ "--debug" ]]; then
			bash fetch.bash $fileout $switchfile "1" "--debug"
		else
			bash fetch.bash $fileout $switchfile "1"
		fi

	fi

	rm fetch.bash
####################		ADD ID's		##########################
	i=1
	total=$(wc -l $switchfile |awk '{print $1}')
	cat $switchfile |while read line
	do
		fasta=$(echo "$line" |awk '{print $1}')
		ti=$(echo "$line" |awk '{print $2}')
		echo "working on $fasta  ($i/$total)"
		headerfasta=$(head -n 1 $fasta)

		if ! [[ "$headerfasta" =~ ti\|[0-9].? ]];then
			sed "s/>/>ti\|$ti\|/g" $fasta > tmp
		else
			cp $fasta tmp
		fi

		case $multiway in
		"0")
			mv tmp new_$fasta
		;;
		"1")
			rm $fasta
			mv tmp $fasta
		;;
		esac

		i=$((i+1))
	done

###################		MERGE FASTAS		############################
	#python merge.py folder_files file_out_name
	echo "* Merging chunk files"
	case $multiway in
	"0")
		rm -f appendheaders.py $fileout $switchfile parsefasta.awk
	;;
	"1")
		makeMergeWork
		$PYTHONBIN merge.py $WORKDIR $multifname.new
		mv $multifname.new new_$multifname
		mv new_$multifname ../.
		cd ..
		rm -rf $multifname""_TMP_FOLDER_DONT_TOUCH $switchfile
	;;
	esac

	echo "* Done :D"
else
	echo "Invalid or Missing Parameters, print --help to see the options"
	exit
fi
