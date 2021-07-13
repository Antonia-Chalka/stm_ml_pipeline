#AMR concantenation
# performed in mobax term
cd /drives/D/PhD/amr_august/concat/amr_out

# copied amr_out to concat as a safe copy

#remove header
#for file in */*.tsv
#do
#  sed -i '1d' $file
#done

#append filename and concantenate
for f in */*.tsv; do
name=${f##*/}
base=${name%????????}
sed "s/$//t $base/" "$f"   # append the number to each line
done > amr_all.tsv

