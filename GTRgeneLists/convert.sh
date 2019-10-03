OIFS="$IFS"
IFS=$'\n'
for f in $(find . -name "*.xlsx")
do
	# echo $f
	# ssconvert -O "separator='     ' format=raw" "$f" "${f%.xlsx}.txt" # Can't change directly to .tsv
	# mv "${f%.xlsx}.txt" "${f%.xlsx}.tsv" # Change the .txt to .tsv
	rm $f # Remove the original Excel worksheets.
done
IFS="$OIFS"
