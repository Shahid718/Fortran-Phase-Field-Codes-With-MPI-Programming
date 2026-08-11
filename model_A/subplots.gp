# file name subplots.gp
set terminal svg size 1200,1000 enhanced font 'Arial,12'
set output 'heatmaps.svg'

# Set up multiplot layout (2 rows, 2 columns)
set multiplot layout 2,2

# Set common color palette
set palette defined (0 "dark-blue", 0.2 "blue", 0.5 "green", 0.8 "yellow", 1 "red")

# Common settings for all plots
# set cbrange [0.49:0.51]  # Adjust based on your data range
set xlabel "Nx"
set ylabel "local Ny"

# Plot 1: ac_rank0.dat
set title "rank0.dat"
plot 'rank0.dat' matrix with image notitle

# Plot 2: ac_rank1.dat
set title "rank1.dat"
plot 'rank1.dat' matrix with image notitle

# Plot 3: ac_rank2.dat
set title "rank2.dat"
plot 'rank2.dat' matrix with image notitle

# Plot 4: ac_rank3.dat
set title "rank3.dat"
plot 'rank3.dat' matrix with image notitle

# Reset multiplot
unset multiplot