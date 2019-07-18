
for i in $(cat temp.txt)
do
  sed -ie "8s/.*/0/" input.dat
  sed -ie "1s/.*/${i}/" input.dat
  ./Monte_Carlo_ISING_1D.exe
  sed -ie "8s/.*/1/" input.dat
  ./Monte_Carlo_ISING_1D.exe
done
