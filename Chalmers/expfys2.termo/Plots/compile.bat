gnuplot timeseries.gpl
epstopdf timeseries.eps
rm timeseries.eps

gnuplot dispersion.gpl
epstopdf dispersion_Al.eps
rm dispersion_Al.eps
#pdflatex dispersion_Al.tex

epstopdf dispersion_Pb.eps
rm dispersion_Pb.eps
#pdflatex dispersion_Pb.tex

gnuplot hyperbolisk_dispersion.gpl
epstopdf hyperbolisk_dispersion.eps
rm hyperbolisk_dispersion.eps
#pdflatex hyperbolisk_dispersion.tex

