Running the mongrail program
----------------------------------------------

Example input files (in example subdirectory):
----------------------------------------------

Chrom filename: c20_m10_r1.chrom
PopA filename: c20_m10_h10_au1_hc0.1.popA
PopB filename: c20_m10_h10_au1_hc0.1.popB
Simulated individual filename: c20_m10_r1_h10_au1_hc0.1.sim


Quick start
---------------------------------------------

git clone
cd mongrail
make; cd example;
../mongrail -p c20_m10_r1.chrom -A c20_m10_h10_au1_hc0.1.popA -B c20_m10_h10_au1_hc0.1.popB -i c20_m10_r1_h10_au1_hc0.1.sim c20_m10_r1_h10_au1_hc0.1.out

