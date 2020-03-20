if [ -e "$1.msh" ]
then
  echo "msh already generated"
else
  ./tools/tetwild/TetWild --save-mid-result 2 $1.stl -e 0.025 --targeted-num-v 1000 --is-laplacian
  mv "$1__mid2.000000.msh" "$1.msh"
fi
