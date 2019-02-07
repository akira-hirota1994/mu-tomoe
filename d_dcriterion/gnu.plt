plot for [i=1:1000] sprintf( "d_txt/meteor_orbit_%d.txt", i ) w l notitle lw 1 lc 7
replot "d_planet/mercury.txt"  w l lw 3 lc 3 title "Mercury"
replot "d_planet/venus.txt"  w l lw 3 lc 5 title "Venus"
replot "d_planet/earth.txt"  w l lw 3 lc 2 title "Earth"
replot "d_planet/mars.txt"  w l lw 3 lc 4 title "Mars"
set key left top
set title "Geminid orbit(Dsh<0.1) num = 819"
set size square
set xl "X[AU]"
set yl "Y[AU]"
set xlabel font "Times-Roman,12"
set ylabel font "Times-Roman,12"
set key font "Times-Roman,10"
set tics font "Times-Roman,10"
set title font "Times-Roman,12"
rep
set terminal pdf
set output "orbit_image.pdf"
rep

