for j in 1 #2 3 4 5
do
for i in /dataH4SM/*
do
	# Set space as the delimiter
IFS='/'

#Read the split words into an array based on space delimiter
read -a strarr <<< "$i"

#Count the total words
#echo "There are ${strarr[3]} words in the text."

	/cytocastCoreNoOsm -i "$i" -o "/outputs/outputH4SM/$j/${strarr[3]}" -t 6 -l 1 -O "/outputs/H4SM/$j/${strarr[3]}";
   #/cytocastCoreNoOsm -i "$i" -o "/app/outputSimplet6steped/$j/${strarr[3]}" -t 6 -l 1 -O "/app/stepLogst6";
done
done
#docker run -it --entrypoint bash --mount type=bind,source="$pwd"/inputFiles/data,target=/app core
