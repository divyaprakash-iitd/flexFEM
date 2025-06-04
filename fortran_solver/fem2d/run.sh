eval "$(conda shell.bash hook)"
conda activate cilia
#python equi_ellip.py
./clean.sh
./newcompile.sh
./ibmc
python plot_images.py
