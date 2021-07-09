folder="output_0/"
#set terminal postfile 
#set title 'Immune deficiency - proliferation\_rate = 2e-3 and death\_rate = 1e-5'
#set logscale y
set ylabel "Number of cell"
set xlabel "Time (hours)"
plot folder.'Cells.dat' u 1:2 w l lw 4 lc 6 title 'lung', folder.'Cells.dat' u 1:3 w l lw 4 lc 5 title 'melanoma', folder.'Cells.dat' u 1:4 w l lw 4 lc 7 title 'CD8', folder.'Cells.dat' u 1:5 w l lw 4 lc 2 title 'macrophage', folder.'Cells.dat' u 1:6 w l lw 4 lc rgb '#FF33CC' title 'dendritic', folder.'Cells.dat' u 1:7 w l lw 4 lc rgb '#FFA500' title 'CD4', folder.'Cells.dat' u 1:8 w l lw 4 lc 8 title 'dead lung', folder.'Cells.dat' u 1:9 w l lw 4 lc 'gray' title 'dead melanoma', folder.'Cells.dat' u 1:10 w l lw 4 lc rgb '#bc8f8f' title 'dead immune'

pause -1 "Hit any key to continue"
#clear

set y2tics
set y2label "TC, TCt, and Tht scale"
set ytics nomirror
set key left top

plot folder.'dm_tc.dat' u 0:1 axis x1y1 w l lw 4 lc 1 title 'DM', folder.'dm_tc.dat' u 0:2 axis x1y2 w l lw 4 lc 2 title 'TC', folder.'dm_tc.dat' u 0:3 axis x1y1 w l lw 4 lc 3 title 'TH1', folder.'dm_tc.dat' u 0:4 axis x1y1 w l lw 4 lc 4 title 'TH2', folder.'dm_tc.dat' u 0:5 axis x1y2 w l lw 4 lc rgb '#FF33CC' title 'TCt', folder.'dm_tc.dat' u 0:6 axis x1y2 w l lw 4 lc 6 title 'Tht'

#plot folder.'dm_tc.dat' u 1:2 w l lw 4 lc 1 title 'DM', folder.'dm_tc.dat' u 1:3 w l lw 4 lc 2 title 'TC', folder.'dm_tc.dat' u 1:4 w l lw 4 lc 3 title 'TH1', folder.'dm_tc.dat' u 1:5 w l lw 4 lc 4 title 'TH2', folder.'dm_tc.dat' u 1:6 w l lw 4 lc rgb '#FF33CC' title 'TCt', folder.'dm_tc.dat' u 1:7 w l lw 4 lc 6 title 'Tht'

pause -1 "Hit any key to continue"
