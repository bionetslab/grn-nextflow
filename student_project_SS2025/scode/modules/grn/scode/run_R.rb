require 'fileutils'

# Input arguments
fdata  = ARGV[0]
ftime  = ARGV[1]
dir    = ARGV[2]
tfnum  = ARGV[3]
pnum   = ARGV[4]
cnum   = ARGV[5]
maxite = ARGV[6]
repnum = ARGV[7].to_i

# Get absolute paths for R scripts
scode_r = File.expand_path(File.join(__dir__, 'SCODE.R'))
avg_r   = File.expand_path(File.join(__dir__, 'averageA.R'))

# Ensure parent output dir exists
FileUtils.mkdir_p(dir)

# Run SCODE multiple times, each in its own subdirectory: results/out_1, out_2, ...
(1..repnum).each do |i|
  subdir = File.join(dir, "out_#{i}")
  FileUtils.mkdir_p(subdir)

  puts "Running SCODE trial #{i}/#{repnum}"
  system("Rscript #{scode_r} #{fdata} #{ftime} #{subdir} #{tfnum} #{pnum} #{cnum} #{maxite}")
end

# Average A matrices
puts "Averaging A matrices"
mean_a_file = File.join(dir, "meanA.txt")
system("Rscript #{avg_r} #{dir} #{repnum} #{mean_a_file}")

# Create ranked edges CSV file
puts "Creating ranked edges file"
ranked_edges_script = File.expand_path(File.join(__dir__, 'create_ranked_edges.R'))
ranked_edges_file = File.join(dir, "ranked_edges.csv")
system("Rscript #{ranked_edges_script} #{mean_a_file} #{ranked_edges_file}")