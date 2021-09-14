

for f in `ls example/input/*.gb`
do
b=`basename $f | sed s/\.gb/\.peg/g`

echo $f $b

python3 gb2peg.py $f example/input/$b

done
