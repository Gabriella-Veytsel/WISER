#!/bin/bash
#SBATCH --job-name=move_rename
#SBATCH --partition=iob_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2gb
#SBATCH --time=1:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@uga.edu


cd /scratch/mandyh/WISER_WW/All_WW_StrainSort

input1='/scratch/mandyh/WISER_WW/WISER_WW_08.17.2023/paired_trim_reads'
input2='/scratch/mandyh/WISER_WW/WISER_WW_11.08.2023/paired_trim_reads'
output='/scratch/mandyh/WISER_WW/All_WW_StrainSort/paired_reads'

cp	$input1/W132_pair_trim_R1.fastq.gz	$output/W1_pair_trim_R1.fastq.gz
cp	$input1/W202_pair_trim_R1.fastq.gz	$output/W2_pair_trim_R1.fastq.gz
cp	$input1/W42_pair_trim_R1.fastq.gz	$output/W3_pair_trim_R1.fastq.gz
cp	$input1/W139_pair_trim_R1.fastq.gz	$output/W4_pair_trim_R1.fastq.gz
cp	$input1/W41_pair_trim_R1.fastq.gz	$output/W5_pair_trim_R1.fastq.gz
cp	$input1/W39_pair_trim_R1.fastq.gz	$output/W6_pair_trim_R1.fastq.gz
cp	$input1/W37_pair_trim_R1.fastq.gz	$output/W7_pair_trim_R1.fastq.gz
cp	$input1/W121_pair_trim_R1.fastq.gz	$output/W8_pair_trim_R1.fastq.gz
cp	$input1/W140_pair_trim_R1.fastq.gz	$output/W9_pair_trim_R1.fastq.gz
cp	$input1/W1_pair_trim_R1.fastq.gz	$output/W10_pair_trim_R1.fastq.gz
cp	$input1/W33_pair_trim_R1.fastq.gz	$output/W11_pair_trim_R1.fastq.gz
cp	$input1/W38_pair_trim_R1.fastq.gz	$output/W12_pair_trim_R1.fastq.gz
cp	$input1/W190_pair_trim_R1.fastq.gz	$output/W13_pair_trim_R1.fastq.gz
cp	$input1/W143_pair_trim_R1.fastq.gz	$output/W14_pair_trim_R1.fastq.gz
cp	$input1/W110_pair_trim_R1.fastq.gz	$output/W15_pair_trim_R1.fastq.gz
cp	$input1/W204_pair_trim_R1.fastq.gz	$output/W16_pair_trim_R1.fastq.gz
cp	$input1/W2_pair_trim_R1.fastq.gz	$output/W17_pair_trim_R1.fastq.gz
cp	$input1/W95_pair_trim_R1.fastq.gz	$output/W18_pair_trim_R1.fastq.gz
cp	$input1/W96_pair_trim_R1.fastq.gz	$output/W19_pair_trim_R1.fastq.gz
cp	$input1/W32_pair_trim_R1.fastq.gz	$output/W20_pair_trim_R1.fastq.gz
cp	$input1/W47_pair_trim_R1.fastq.gz	$output/W21_pair_trim_R1.fastq.gz
cp	$input1/W178_pair_trim_R1.fastq.gz	$output/W22_pair_trim_R1.fastq.gz
cp	$input1/W169_pair_trim_R1.fastq.gz	$output/W23_pair_trim_R1.fastq.gz
cp	$input1/W205_pair_trim_R1.fastq.gz	$output/W24_pair_trim_R1.fastq.gz
cp	$input1/W14_pair_trim_R1.fastq.gz	$output/W25_pair_trim_R1.fastq.gz
cp	$input1/W152_pair_trim_R1.fastq.gz	$output/W26_pair_trim_R1.fastq.gz
cp	$input1/W22_pair_trim_R1.fastq.gz	$output/W27_pair_trim_R1.fastq.gz
cp	$input1/W191_pair_trim_R1.fastq.gz	$output/W28_pair_trim_R1.fastq.gz
cp	$input1/W136_pair_trim_R1.fastq.gz	$output/W29_pair_trim_R1.fastq.gz
cp	$input1/W45_pair_trim_R1.fastq.gz	$output/W30_pair_trim_R1.fastq.gz
cp	$input1/W94_pair_trim_R1.fastq.gz	$output/W31_pair_trim_R1.fastq.gz
cp	$input1/W118_pair_trim_R1.fastq.gz	$output/W32_pair_trim_R1.fastq.gz
cp	$input1/W20_pair_trim_R1.fastq.gz	$output/W33_pair_trim_R1.fastq.gz
cp	$input1/W12_pair_trim_R1.fastq.gz	$output/W34_pair_trim_R1.fastq.gz
cp	$input1/W134_pair_trim_R1.fastq.gz	$output/W35_pair_trim_R1.fastq.gz
cp	$input1/W201_pair_trim_R1.fastq.gz	$output/W36_pair_trim_R1.fastq.gz
cp	$input1/W124_pair_trim_R1.fastq.gz	$output/W37_pair_trim_R1.fastq.gz
cp	$input1/W97_pair_trim_R1.fastq.gz	$output/W38_pair_trim_R1.fastq.gz
cp	$input1/W174_pair_trim_R1.fastq.gz	$output/W39_pair_trim_R1.fastq.gz
cp	$input1/W25_pair_trim_R1.fastq.gz	$output/W40_pair_trim_R1.fastq.gz
cp	$input1/W4_pair_trim_R1.fastq.gz	$output/W41_pair_trim_R1.fastq.gz
cp	$input1/W188_pair_trim_R1.fastq.gz	$output/W42_pair_trim_R1.fastq.gz
cp	$input1/W171_pair_trim_R1.fastq.gz	$output/W43_pair_trim_R1.fastq.gz
cp	$input1/W102_pair_trim_R1.fastq.gz	$output/W44_pair_trim_R1.fastq.gz
cp	$input1/W172_pair_trim_R1.fastq.gz	$output/W45_pair_trim_R1.fastq.gz
cp	$input1/W28_pair_trim_R1.fastq.gz	$output/W46_pair_trim_R1.fastq.gz
cp	$input1/W34_pair_trim_R1.fastq.gz	$output/W47_pair_trim_R1.fastq.gz
cp	$input1/W23_pair_trim_R1.fastq.gz	$output/W48_pair_trim_R1.fastq.gz
cp	$input1/W176_pair_trim_R1.fastq.gz	$output/W49_pair_trim_R1.fastq.gz
cp	$input1/W11_pair_trim_R1.fastq.gz	$output/W50_pair_trim_R1.fastq.gz
cp	$input1/W119_pair_trim_R1.fastq.gz	$output/W51_pair_trim_R1.fastq.gz
cp	$input1/W105_pair_trim_R1.fastq.gz	$output/W52_pair_trim_R1.fastq.gz
cp	$input1/W17_pair_trim_R1.fastq.gz	$output/W53_pair_trim_R1.fastq.gz
cp	$input1/W18_pair_trim_R1.fastq.gz	$output/W54_pair_trim_R1.fastq.gz
cp	$input1/W99_pair_trim_R1.fastq.gz	$output/W55_pair_trim_R1.fastq.gz
cp	$input1/W13_pair_trim_R1.fastq.gz	$output/W56_pair_trim_R1.fastq.gz
cp	$input1/W151_pair_trim_R1.fastq.gz	$output/W57_pair_trim_R1.fastq.gz
cp	$input1/W141_pair_trim_R1.fastq.gz	$output/W58_pair_trim_R1.fastq.gz
cp	$input1/W144_pair_trim_R1.fastq.gz	$output/W59_pair_trim_R1.fastq.gz
cp	$input1/W6_pair_trim_R1.fastq.gz	$output/W60_pair_trim_R1.fastq.gz
cp	$input1/W145_pair_trim_R1.fastq.gz	$output/W61_pair_trim_R1.fastq.gz
cp	$input1/W111_pair_trim_R1.fastq.gz	$output/W62_pair_trim_R1.fastq.gz
cp	$input1/W16_pair_trim_R1.fastq.gz	$output/W63_pair_trim_R1.fastq.gz
cp	$input1/W21_pair_trim_R1.fastq.gz	$output/W64_pair_trim_R1.fastq.gz
cp	$input1/W49_pair_trim_R1.fastq.gz	$output/W65_pair_trim_R1.fastq.gz
cp	$input1/W100_pair_trim_R1.fastq.gz	$output/W66_pair_trim_R1.fastq.gz
cp	$input1/W138_pair_trim_R1.fastq.gz	$output/W67_pair_trim_R1.fastq.gz
cp	$input1/W30_pair_trim_R1.fastq.gz	$output/W68_pair_trim_R1.fastq.gz
cp	$input1/W7_pair_trim_R1.fastq.gz	$output/W69_pair_trim_R1.fastq.gz
cp	$input1/W133_pair_trim_R1.fastq.gz	$output/W70_pair_trim_R1.fastq.gz
cp	$input1/W48_pair_trim_R1.fastq.gz	$output/W71_pair_trim_R1.fastq.gz
cp	$input1/W193_pair_trim_R1.fastq.gz	$output/W72_pair_trim_R1.fastq.gz
cp	$input1/W5_pair_trim_R1.fastq.gz	$output/W73_pair_trim_R1.fastq.gz
cp	$input1/W26_pair_trim_R1.fastq.gz	$output/W74_pair_trim_R1.fastq.gz
cp	$input1/W19_pair_trim_R1.fastq.gz	$output/W75_pair_trim_R1.fastq.gz
cp	$input1/W106_pair_trim_R1.fastq.gz	$output/W76_pair_trim_R1.fastq.gz
cp	$input1/W3_pair_trim_R1.fastq.gz	$output/W77_pair_trim_R1.fastq.gz
cp	$input1/W149_pair_trim_R1.fastq.gz	$output/W78_pair_trim_R1.fastq.gz
cp	$input1/W93_pair_trim_R1.fastq.gz	$output/W79_pair_trim_R1.fastq.gz
cp	$input1/W122_pair_trim_R1.fastq.gz	$output/W80_pair_trim_R1.fastq.gz
cp	$input1/W194_pair_trim_R1.fastq.gz	$output/W81_pair_trim_R1.fastq.gz
cp	$input1/W173_pair_trim_R1.fastq.gz	$output/W82_pair_trim_R1.fastq.gz
cp	$input1/W36_pair_trim_R1.fastq.gz	$output/W83_pair_trim_R1.fastq.gz
cp	$input1/W179_pair_trim_R1.fastq.gz	$output/W84_pair_trim_R1.fastq.gz
cp	$input1/W164_pair_trim_R1.fastq.gz	$output/W85_pair_trim_R1.fastq.gz
cp	$input1/W44_pair_trim_R1.fastq.gz	$output/W86_pair_trim_R1.fastq.gz
cp	$input1/W185_pair_trim_R1.fastq.gz	$output/W87_pair_trim_R1.fastq.gz
cp	$input1/W24_pair_trim_R1.fastq.gz	$output/W88_pair_trim_R1.fastq.gz
cp	$input1/W117_pair_trim_R1.fastq.gz	$output/W89_pair_trim_R1.fastq.gz
cp	$input1/W197_pair_trim_R1.fastq.gz	$output/W90_pair_trim_R1.fastq.gz
cp	$input1/W147_pair_trim_R1.fastq.gz	$output/W91_pair_trim_R1.fastq.gz
cp	$input1/W31_pair_trim_R1.fastq.gz	$output/W92_pair_trim_R1.fastq.gz
cp	$input1/W10_pair_trim_R1.fastq.gz	$output/W93_pair_trim_R1.fastq.gz
cp	$input1/W123_pair_trim_R1.fastq.gz	$output/W94_pair_trim_R1.fastq.gz
cp	$input1/W9_pair_trim_R1.fastq.gz	$output/W95_pair_trim_R1.fastq.gz
cp	$input1/W15_pair_trim_R1.fastq.gz	$output/W96_pair_trim_R1.fastq.gz
cp	$input1/W40_pair_trim_R1.fastq.gz	$output/W97_pair_trim_R1.fastq.gz
cp	$input1/W35_pair_trim_R1.fastq.gz	$output/W98_pair_trim_R1.fastq.gz
cp	$input1/W101_pair_trim_R1.fastq.gz	$output/W99_pair_trim_R1.fastq.gz
cp	$input1/W142_pair_trim_R1.fastq.gz	$output/W100_pair_trim_R1.fastq.gz
cp	$input1/W130_pair_trim_R1.fastq.gz	$output/W101_pair_trim_R1.fastq.gz
cp	$input1/W29_pair_trim_R1.fastq.gz	$output/W102_pair_trim_R1.fastq.gz
cp	$input1/W189_pair_trim_R1.fastq.gz	$output/W103_pair_trim_R1.fastq.gz
cp	$input1/W163_pair_trim_R1.fastq.gz	$output/W104_pair_trim_R1.fastq.gz
cp	$input1/W108_pair_trim_R1.fastq.gz	$output/W105_pair_trim_R1.fastq.gz
cp	$input1/W27_pair_trim_R1.fastq.gz	$output/W106_pair_trim_R1.fastq.gz
cp	$input1/W74_pair_trim_R1.fastq.gz	$output/W107_pair_trim_R1.fastq.gz
cp	$input1/W196_pair_trim_R1.fastq.gz	$output/W108_pair_trim_R1.fastq.gz
cp	$input1/W146_pair_trim_R1.fastq.gz	$output/W109_pair_trim_R1.fastq.gz
cp	$input1/W162_pair_trim_R1.fastq.gz	$output/W110_pair_trim_R1.fastq.gz
cp	$input1/W131_pair_trim_R1.fastq.gz	$output/W111_pair_trim_R1.fastq.gz
cp	$input1/W63_pair_trim_R1.fastq.gz	$output/W112_pair_trim_R1.fastq.gz
cp	$input1/W180_pair_trim_R1.fastq.gz	$output/W113_pair_trim_R1.fastq.gz
cp	$input1/W150_pair_trim_R1.fastq.gz	$output/W114_pair_trim_R1.fastq.gz
cp	$input1/W203_pair_trim_R1.fastq.gz	$output/W115_pair_trim_R1.fastq.gz
cp	$input1/W137_pair_trim_R1.fastq.gz	$output/W116_pair_trim_R1.fastq.gz
cp	$input1/W128_pair_trim_R1.fastq.gz	$output/W117_pair_trim_R1.fastq.gz
cp	$input1/W165_pair_trim_R1.fastq.gz	$output/W118_pair_trim_R1.fastq.gz
cp	$input1/W103_pair_trim_R1.fastq.gz	$output/W119_pair_trim_R1.fastq.gz
cp	$input1/W195_pair_trim_R1.fastq.gz	$output/W120_pair_trim_R1.fastq.gz
cp	$input1/W120_pair_trim_R1.fastq.gz	$output/W121_pair_trim_R1.fastq.gz
cp	$input1/W127_pair_trim_R1.fastq.gz	$output/W122_pair_trim_R1.fastq.gz
cp	$input1/W129_pair_trim_R1.fastq.gz	$output/W123_pair_trim_R1.fastq.gz
cp	$input1/W192_pair_trim_R1.fastq.gz	$output/W124_pair_trim_R1.fastq.gz
cp	$input1/W166_pair_trim_R1.fastq.gz	$output/W125_pair_trim_R1.fastq.gz
cp	$input1/W159_pair_trim_R1.fastq.gz	$output/W126_pair_trim_R1.fastq.gz
cp	$input1/W177_pair_trim_R1.fastq.gz	$output/W127_pair_trim_R1.fastq.gz
cp	$input1/W126_pair_trim_R1.fastq.gz	$output/W128_pair_trim_R1.fastq.gz
cp	$input1/W183_pair_trim_R1.fastq.gz	$output/W129_pair_trim_R1.fastq.gz
cp	$input1/W98_pair_trim_R1.fastq.gz	$output/W130_pair_trim_R1.fastq.gz
cp	$input1/W132_pair_trim_R2.fastq.gz	$output/W1_pair_trim_R2.fastq.gz
cp	$input1/W202_pair_trim_R2.fastq.gz	$output/W2_pair_trim_R2.fastq.gz
cp	$input1/W42_pair_trim_R2.fastq.gz	$output/W3_pair_trim_R2.fastq.gz
cp	$input1/W139_pair_trim_R2.fastq.gz	$output/W4_pair_trim_R2.fastq.gz
cp	$input1/W41_pair_trim_R2.fastq.gz	$output/W5_pair_trim_R2.fastq.gz
cp	$input1/W39_pair_trim_R2.fastq.gz	$output/W6_pair_trim_R2.fastq.gz
cp	$input1/W37_pair_trim_R2.fastq.gz	$output/W7_pair_trim_R2.fastq.gz
cp	$input1/W121_pair_trim_R2.fastq.gz	$output/W8_pair_trim_R2.fastq.gz
cp	$input1/W140_pair_trim_R2.fastq.gz	$output/W9_pair_trim_R2.fastq.gz
cp	$input1/W1_pair_trim_R2.fastq.gz	$output/W10_pair_trim_R2.fastq.gz
cp	$input1/W33_pair_trim_R2.fastq.gz	$output/W11_pair_trim_R2.fastq.gz
cp	$input1/W38_pair_trim_R2.fastq.gz	$output/W12_pair_trim_R2.fastq.gz
cp	$input1/W190_pair_trim_R2.fastq.gz	$output/W13_pair_trim_R2.fastq.gz
cp	$input1/W143_pair_trim_R2.fastq.gz	$output/W14_pair_trim_R2.fastq.gz
cp	$input1/W110_pair_trim_R2.fastq.gz	$output/W15_pair_trim_R2.fastq.gz
cp	$input1/W204_pair_trim_R2.fastq.gz	$output/W16_pair_trim_R2.fastq.gz
cp	$input1/W2_pair_trim_R2.fastq.gz	$output/W17_pair_trim_R2.fastq.gz
cp	$input1/W95_pair_trim_R2.fastq.gz	$output/W18_pair_trim_R2.fastq.gz
cp	$input1/W96_pair_trim_R2.fastq.gz	$output/W19_pair_trim_R2.fastq.gz
cp	$input1/W32_pair_trim_R2.fastq.gz	$output/W20_pair_trim_R2.fastq.gz
cp	$input1/W47_pair_trim_R2.fastq.gz	$output/W21_pair_trim_R2.fastq.gz
cp	$input1/W178_pair_trim_R2.fastq.gz	$output/W22_pair_trim_R2.fastq.gz
cp	$input1/W169_pair_trim_R2.fastq.gz	$output/W23_pair_trim_R2.fastq.gz
cp	$input1/W205_pair_trim_R2.fastq.gz	$output/W24_pair_trim_R2.fastq.gz
cp	$input1/W14_pair_trim_R2.fastq.gz	$output/W25_pair_trim_R2.fastq.gz
cp	$input1/W152_pair_trim_R2.fastq.gz	$output/W26_pair_trim_R2.fastq.gz
cp	$input1/W22_pair_trim_R2.fastq.gz	$output/W27_pair_trim_R2.fastq.gz
cp	$input1/W191_pair_trim_R2.fastq.gz	$output/W28_pair_trim_R2.fastq.gz
cp	$input1/W136_pair_trim_R2.fastq.gz	$output/W29_pair_trim_R2.fastq.gz
cp	$input1/W45_pair_trim_R2.fastq.gz	$output/W30_pair_trim_R2.fastq.gz
cp	$input1/W94_pair_trim_R2.fastq.gz	$output/W31_pair_trim_R2.fastq.gz
cp	$input1/W118_pair_trim_R2.fastq.gz	$output/W32_pair_trim_R2.fastq.gz
cp	$input1/W20_pair_trim_R2.fastq.gz	$output/W33_pair_trim_R2.fastq.gz
cp	$input1/W12_pair_trim_R2.fastq.gz	$output/W34_pair_trim_R2.fastq.gz
cp	$input1/W134_pair_trim_R2.fastq.gz	$output/W35_pair_trim_R2.fastq.gz
cp	$input1/W201_pair_trim_R2.fastq.gz	$output/W36_pair_trim_R2.fastq.gz
cp	$input1/W124_pair_trim_R2.fastq.gz	$output/W37_pair_trim_R2.fastq.gz
cp	$input1/W97_pair_trim_R2.fastq.gz	$output/W38_pair_trim_R2.fastq.gz
cp	$input1/W174_pair_trim_R2.fastq.gz	$output/W39_pair_trim_R2.fastq.gz
cp	$input1/W25_pair_trim_R2.fastq.gz	$output/W40_pair_trim_R2.fastq.gz
cp	$input1/W4_pair_trim_R2.fastq.gz	$output/W41_pair_trim_R2.fastq.gz
cp	$input1/W188_pair_trim_R2.fastq.gz	$output/W42_pair_trim_R2.fastq.gz
cp	$input1/W171_pair_trim_R2.fastq.gz	$output/W43_pair_trim_R2.fastq.gz
cp	$input1/W102_pair_trim_R2.fastq.gz	$output/W44_pair_trim_R2.fastq.gz
cp	$input1/W172_pair_trim_R2.fastq.gz	$output/W45_pair_trim_R2.fastq.gz
cp	$input1/W28_pair_trim_R2.fastq.gz	$output/W46_pair_trim_R2.fastq.gz
cp	$input1/W34_pair_trim_R2.fastq.gz	$output/W47_pair_trim_R2.fastq.gz
cp	$input1/W23_pair_trim_R2.fastq.gz	$output/W48_pair_trim_R2.fastq.gz
cp	$input1/W176_pair_trim_R2.fastq.gz	$output/W49_pair_trim_R2.fastq.gz
cp	$input1/W11_pair_trim_R2.fastq.gz	$output/W50_pair_trim_R2.fastq.gz
cp	$input1/W119_pair_trim_R2.fastq.gz	$output/W51_pair_trim_R2.fastq.gz
cp	$input1/W105_pair_trim_R2.fastq.gz	$output/W52_pair_trim_R2.fastq.gz
cp	$input1/W17_pair_trim_R2.fastq.gz	$output/W53_pair_trim_R2.fastq.gz
cp	$input1/W18_pair_trim_R2.fastq.gz	$output/W54_pair_trim_R2.fastq.gz
cp	$input1/W99_pair_trim_R2.fastq.gz	$output/W55_pair_trim_R2.fastq.gz
cp	$input1/W13_pair_trim_R2.fastq.gz	$output/W56_pair_trim_R2.fastq.gz
cp	$input1/W151_pair_trim_R2.fastq.gz	$output/W57_pair_trim_R2.fastq.gz
cp	$input1/W141_pair_trim_R2.fastq.gz	$output/W58_pair_trim_R2.fastq.gz
cp	$input1/W144_pair_trim_R2.fastq.gz	$output/W59_pair_trim_R2.fastq.gz
cp	$input1/W6_pair_trim_R2.fastq.gz	$output/W60_pair_trim_R2.fastq.gz
cp	$input1/W145_pair_trim_R2.fastq.gz	$output/W61_pair_trim_R2.fastq.gz
cp	$input1/W111_pair_trim_R2.fastq.gz	$output/W62_pair_trim_R2.fastq.gz
cp	$input1/W16_pair_trim_R2.fastq.gz	$output/W63_pair_trim_R2.fastq.gz
cp	$input1/W21_pair_trim_R2.fastq.gz	$output/W64_pair_trim_R2.fastq.gz
cp	$input1/W49_pair_trim_R2.fastq.gz	$output/W65_pair_trim_R2.fastq.gz
cp	$input1/W100_pair_trim_R2.fastq.gz	$output/W66_pair_trim_R2.fastq.gz
cp	$input1/W138_pair_trim_R2.fastq.gz	$output/W67_pair_trim_R2.fastq.gz
cp	$input1/W30_pair_trim_R2.fastq.gz	$output/W68_pair_trim_R2.fastq.gz
cp	$input1/W7_pair_trim_R2.fastq.gz	$output/W69_pair_trim_R2.fastq.gz
cp	$input1/W133_pair_trim_R2.fastq.gz	$output/W70_pair_trim_R2.fastq.gz
cp	$input1/W48_pair_trim_R2.fastq.gz	$output/W71_pair_trim_R2.fastq.gz
cp	$input1/W193_pair_trim_R2.fastq.gz	$output/W72_pair_trim_R2.fastq.gz
cp	$input1/W5_pair_trim_R2.fastq.gz	$output/W73_pair_trim_R2.fastq.gz
cp	$input1/W26_pair_trim_R2.fastq.gz	$output/W74_pair_trim_R2.fastq.gz
cp	$input1/W19_pair_trim_R2.fastq.gz	$output/W75_pair_trim_R2.fastq.gz
cp	$input1/W106_pair_trim_R2.fastq.gz	$output/W76_pair_trim_R2.fastq.gz
cp	$input1/W3_pair_trim_R2.fastq.gz	$output/W77_pair_trim_R2.fastq.gz
cp	$input1/W149_pair_trim_R2.fastq.gz	$output/W78_pair_trim_R2.fastq.gz
cp	$input1/W93_pair_trim_R2.fastq.gz	$output/W79_pair_trim_R2.fastq.gz
cp	$input1/W122_pair_trim_R2.fastq.gz	$output/W80_pair_trim_R2.fastq.gz
cp	$input1/W194_pair_trim_R2.fastq.gz	$output/W81_pair_trim_R2.fastq.gz
cp	$input1/W173_pair_trim_R2.fastq.gz	$output/W82_pair_trim_R2.fastq.gz
cp	$input1/W36_pair_trim_R2.fastq.gz	$output/W83_pair_trim_R2.fastq.gz
cp	$input1/W179_pair_trim_R2.fastq.gz	$output/W84_pair_trim_R2.fastq.gz
cp	$input1/W164_pair_trim_R2.fastq.gz	$output/W85_pair_trim_R2.fastq.gz
cp	$input1/W44_pair_trim_R2.fastq.gz	$output/W86_pair_trim_R2.fastq.gz
cp	$input1/W185_pair_trim_R2.fastq.gz	$output/W87_pair_trim_R2.fastq.gz
cp	$input1/W24_pair_trim_R2.fastq.gz	$output/W88_pair_trim_R2.fastq.gz
cp	$input1/W117_pair_trim_R2.fastq.gz	$output/W89_pair_trim_R2.fastq.gz
cp	$input1/W197_pair_trim_R2.fastq.gz	$output/W90_pair_trim_R2.fastq.gz
cp	$input1/W147_pair_trim_R2.fastq.gz	$output/W91_pair_trim_R2.fastq.gz
cp	$input1/W31_pair_trim_R2.fastq.gz	$output/W92_pair_trim_R2.fastq.gz
cp	$input1/W10_pair_trim_R2.fastq.gz	$output/W93_pair_trim_R2.fastq.gz
cp	$input1/W123_pair_trim_R2.fastq.gz	$output/W94_pair_trim_R2.fastq.gz
cp	$input1/W9_pair_trim_R2.fastq.gz	$output/W95_pair_trim_R2.fastq.gz
cp	$input1/W15_pair_trim_R2.fastq.gz	$output/W96_pair_trim_R2.fastq.gz
cp	$input1/W40_pair_trim_R2.fastq.gz	$output/W97_pair_trim_R2.fastq.gz
cp	$input1/W35_pair_trim_R2.fastq.gz	$output/W98_pair_trim_R2.fastq.gz
cp	$input1/W101_pair_trim_R2.fastq.gz	$output/W99_pair_trim_R2.fastq.gz
cp	$input1/W142_pair_trim_R2.fastq.gz	$output/W100_pair_trim_R2.fastq.gz
cp	$input1/W130_pair_trim_R2.fastq.gz	$output/W101_pair_trim_R2.fastq.gz
cp	$input1/W29_pair_trim_R2.fastq.gz	$output/W102_pair_trim_R2.fastq.gz
cp	$input1/W189_pair_trim_R2.fastq.gz	$output/W103_pair_trim_R2.fastq.gz
cp	$input1/W163_pair_trim_R2.fastq.gz	$output/W104_pair_trim_R2.fastq.gz
cp	$input1/W108_pair_trim_R2.fastq.gz	$output/W105_pair_trim_R2.fastq.gz
cp	$input1/W27_pair_trim_R2.fastq.gz	$output/W106_pair_trim_R2.fastq.gz
cp	$input1/W74_pair_trim_R2.fastq.gz	$output/W107_pair_trim_R2.fastq.gz
cp	$input1/W196_pair_trim_R2.fastq.gz	$output/W108_pair_trim_R2.fastq.gz
cp	$input1/W146_pair_trim_R2.fastq.gz	$output/W109_pair_trim_R2.fastq.gz
cp	$input1/W162_pair_trim_R2.fastq.gz	$output/W110_pair_trim_R2.fastq.gz
cp	$input1/W131_pair_trim_R2.fastq.gz	$output/W111_pair_trim_R2.fastq.gz
cp	$input1/W63_pair_trim_R2.fastq.gz	$output/W112_pair_trim_R2.fastq.gz
cp	$input1/W180_pair_trim_R2.fastq.gz	$output/W113_pair_trim_R2.fastq.gz
cp	$input1/W150_pair_trim_R2.fastq.gz	$output/W114_pair_trim_R2.fastq.gz
cp	$input1/W203_pair_trim_R2.fastq.gz	$output/W115_pair_trim_R2.fastq.gz
cp	$input1/W137_pair_trim_R2.fastq.gz	$output/W116_pair_trim_R2.fastq.gz
cp	$input1/W128_pair_trim_R2.fastq.gz	$output/W117_pair_trim_R2.fastq.gz
cp	$input1/W165_pair_trim_R2.fastq.gz	$output/W118_pair_trim_R2.fastq.gz
cp	$input1/W103_pair_trim_R2.fastq.gz	$output/W119_pair_trim_R2.fastq.gz
cp	$input1/W195_pair_trim_R2.fastq.gz	$output/W120_pair_trim_R2.fastq.gz
cp	$input1/W120_pair_trim_R2.fastq.gz	$output/W121_pair_trim_R2.fastq.gz
cp	$input1/W127_pair_trim_R2.fastq.gz	$output/W122_pair_trim_R2.fastq.gz
cp	$input1/W129_pair_trim_R2.fastq.gz	$output/W123_pair_trim_R2.fastq.gz
cp	$input1/W192_pair_trim_R2.fastq.gz	$output/W124_pair_trim_R2.fastq.gz
cp	$input1/W166_pair_trim_R2.fastq.gz	$output/W125_pair_trim_R2.fastq.gz
cp	$input1/W159_pair_trim_R2.fastq.gz	$output/W126_pair_trim_R2.fastq.gz
cp	$input1/W177_pair_trim_R2.fastq.gz	$output/W127_pair_trim_R2.fastq.gz
cp	$input1/W126_pair_trim_R2.fastq.gz	$output/W128_pair_trim_R2.fastq.gz
cp	$input1/W183_pair_trim_R2.fastq.gz	$output/W129_pair_trim_R2.fastq.gz
cp	$input1/W98_pair_trim_R2.fastq.gz	$output/W130_pair_trim_R2.fastq.gz
cp	$input2/COLL_73_S39_pair_trim_R1.fastq.gz	$output/W131_pair_trim_R1.fastq.gz
cp	$input2/COLL_238_S34_pair_trim_R1.fastq.gz	$output/W132_pair_trim_R1.fastq.gz
cp	$input2/COLL_240_S35_pair_trim_R1.fastq.gz	$output/W133_pair_trim_R1.fastq.gz
cp	$input2/COLL_217_S21_pair_trim_R1.fastq.gz	$output/W134_pair_trim_R1.fastq.gz
cp	$input2/COLL_100_S37_pair_trim_R1.fastq.gz	$output/W135_pair_trim_R1.fastq.gz
cp	$input2/NO_230_S47_pair_trim_R1.fastq.gz	$output/W136_pair_trim_R1.fastq.gz
cp	$input2/COLL_216_S20_pair_trim_R1.fastq.gz	$output/W137_pair_trim_R1.fastq.gz
cp	$input2/MI_228_S46_pair_trim_R1.fastq.gz	$output/W138_pair_trim_R1.fastq.gz
cp	$input2/MI_122_S28_pair_trim_R1.fastq.gz	$output/W139_pair_trim_R1.fastq.gz
cp	$input2/CC_132_S31_pair_trim_R1.fastq.gz	$output/W140_pair_trim_R1.fastq.gz
cp	$input2/COLL_170_S3_pair_trim_R1.fastq.gz	$output/W141_pair_trim_R1.fastq.gz
cp	$input2/MI_123_S29_pair_trim_R1.fastq.gz	$output/W142_pair_trim_R1.fastq.gz
cp	$input2/COLL_209_S36_pair_trim_R1.fastq.gz	$output/W143_pair_trim_R1.fastq.gz
cp	$input2/COLL_118_S27_pair_trim_R1.fastq.gz	$output/W144_pair_trim_R1.fastq.gz
cp	$input2/MI_233_S51_pair_trim_R1.fastq.gz	$output/W145_pair_trim_R1.fastq.gz
cp	$input2/MI_142_S32_pair_trim_R1.fastq.gz	$output/W146_pair_trim_R1.fastq.gz
cp	$input2/COLL_184_S11_pair_trim_R1.fastq.gz	$output/W147_pair_trim_R1.fastq.gz
cp	$input2/POS_CTRL_6_S40_pair_trim_R1.fastq.gz	$output/W148_pair_trim_R1.fastq.gz
cp	$input2/NO_126_S30_pair_trim_R1.fastq.gz	$output/W149_pair_trim_R1.fastq.gz
cp	$input2/COLL_214_S19_pair_trim_R1.fastq.gz	$output/W150_pair_trim_R1.fastq.gz
cp	$input2/NEG_CTRL_6_S24_pair_trim_R1.fastq.gz	$output/W151_pair_trim_R1.fastq.gz
cp	$input2/NO_171_S4_pair_trim_R1.fastq.gz	$output/W152_pair_trim_R1.fastq.gz
cp	$input2/COLL_210_S18_pair_trim_R1.fastq.gz	$output/W153_pair_trim_R1.fastq.gz
cp	$input2/COLL_220_S25_pair_trim_R1.fastq.gz	$output/W154_pair_trim_R1.fastq.gz
cp	$input2/POS_CTRL_7_S55_pair_trim_R1.fastq.gz	$output/W155_pair_trim_R1.fastq.gz
cp	$input2/COLL_232_S50_pair_trim_R1.fastq.gz	$output/W156_pair_trim_R1.fastq.gz
cp	$input2/COLL_200_S13_pair_trim_R1.fastq.gz	$output/W157_pair_trim_R1.fastq.gz
cp	$input2/COLL_97_S38_pair_trim_R1.fastq.gz	$output/W158_pair_trim_R1.fastq.gz
cp	$input2/COLL_219_S23_pair_trim_R1.fastq.gz	$output/W159_pair_trim_R1.fastq.gz
cp	$input2/COLL_176_S8_pair_trim_R1.fastq.gz	$output/W160_pair_trim_R1.fastq.gz
cp	$input2/COLL_208_S17_pair_trim_R1.fastq.gz	$output/W161_pair_trim_R1.fastq.gz
cp	$input2/COLL_231_S49_pair_trim_R1.fastq.gz	$output/W162_pair_trim_R1.fastq.gz
cp	$input2/COLL_177_S9_pair_trim_R1.fastq.gz	$output/W163_pair_trim_R1.fastq.gz
cp	$input2/COLL_201_S14_pair_trim_R1.fastq.gz	$output/W164_pair_trim_R1.fastq.gz
cp	$input2/COLL_234_S52_pair_trim_R1.fastq.gz	$output/W165_pair_trim_R1.fastq.gz
cp	$input2/NO_147_S33_pair_trim_R1.fastq.gz	$output/W166_pair_trim_R1.fastq.gz
cp	$input2/COLL_205_S16_pair_trim_R1.fastq.gz	$output/W167_pair_trim_R1.fastq.gz
cp	$input2/COLL_179_S10_pair_trim_R1.fastq.gz	$output/W168_pair_trim_R1.fastq.gz
cp	$input2/COLL_174_S6_pair_trim_R1.fastq.gz	$output/W169_pair_trim_R1.fastq.gz
cp	$input2/COLL_160_S1_pair_trim_R1.fastq.gz	$output/W170_pair_trim_R1.fastq.gz
cp	$input2/COLL_173_S5_pair_trim_R1.fastq.gz	$output/W171_pair_trim_R1.fastq.gz
cp	$input2/COLL_226_S44_pair_trim_R1.fastq.gz	$output/W172_pair_trim_R1.fastq.gz
cp	$input2/COLL_57_S41_pair_trim_R1.fastq.gz	$output/W173_pair_trim_R1.fastq.gz
cp	$input2/COLL_218_S22_pair_trim_R1.fastq.gz	$output/W174_pair_trim_R1.fastq.gz
cp	$input2/NEG_CTRL_7_S48_pair_trim_R1.fastq.gz	$output/W175_pair_trim_R1.fastq.gz
cp	$input2/COLL_213_S43_pair_trim_R1.fastq.gz	$output/W176_pair_trim_R1.fastq.gz
cp	$input2/COLL_236_S54_pair_trim_R1.fastq.gz	$output/W177_pair_trim_R1.fastq.gz
cp	$input2/COLL_118_S42_pair_trim_R1.fastq.gz	$output/W178_pair_trim_R1.fastq.gz
cp	$input2/COLL_73_S39_pair_trim_R2.fastq.gz	$output/W131_pair_trim_R2.fastq.gz
cp	$input2/COLL_238_S34_pair_trim_R2.fastq.gz	$output/W132_pair_trim_R2.fastq.gz
cp	$input2/COLL_240_S35_pair_trim_R2.fastq.gz	$output/W133_pair_trim_R2.fastq.gz
cp	$input2/COLL_217_S21_pair_trim_R2.fastq.gz	$output/W134_pair_trim_R2.fastq.gz
cp	$input2/COLL_100_S37_pair_trim_R2.fastq.gz	$output/W135_pair_trim_R2.fastq.gz
cp	$input2/NO_230_S47_pair_trim_R2.fastq.gz	$output/W136_pair_trim_R2.fastq.gz
cp	$input2/COLL_216_S20_pair_trim_R2.fastq.gz	$output/W137_pair_trim_R2.fastq.gz
cp	$input2/MI_228_S46_pair_trim_R2.fastq.gz	$output/W138_pair_trim_R2.fastq.gz
cp	$input2/MI_122_S28_pair_trim_R2.fastq.gz	$output/W139_pair_trim_R2.fastq.gz
cp	$input2/CC_132_S31_pair_trim_R2.fastq.gz	$output/W140_pair_trim_R2.fastq.gz
cp	$input2/COLL_170_S3_pair_trim_R2.fastq.gz	$output/W141_pair_trim_R2.fastq.gz
cp	$input2/MI_123_S29_pair_trim_R2.fastq.gz	$output/W142_pair_trim_R2.fastq.gz
cp	$input2/COLL_209_S36_pair_trim_R2.fastq.gz	$output/W143_pair_trim_R2.fastq.gz
cp	$input2/COLL_118_S27_pair_trim_R2.fastq.gz	$output/W144_pair_trim_R2.fastq.gz
cp	$input2/MI_233_S51_pair_trim_R2.fastq.gz	$output/W145_pair_trim_R2.fastq.gz
cp	$input2/MI_142_S32_pair_trim_R2.fastq.gz	$output/W146_pair_trim_R2.fastq.gz
cp	$input2/COLL_184_S11_pair_trim_R2.fastq.gz	$output/W147_pair_trim_R2.fastq.gz
cp	$input2/POS_CTRL_6_S40_pair_trim_R2.fastq.gz	$output/W148_pair_trim_R2.fastq.gz
cp	$input2/NO_126_S30_pair_trim_R2.fastq.gz	$output/W149_pair_trim_R2.fastq.gz
cp	$input2/COLL_214_S19_pair_trim_R2.fastq.gz	$output/W150_pair_trim_R2.fastq.gz
cp	$input2/NEG_CTRL_6_S24_pair_trim_R2.fastq.gz	$output/W151_pair_trim_R2.fastq.gz
cp	$input2/NO_171_S4_pair_trim_R2.fastq.gz	$output/W152_pair_trim_R2.fastq.gz
cp	$input2/COLL_210_S18_pair_trim_R2.fastq.gz	$output/W153_pair_trim_R2.fastq.gz
cp	$input2/COLL_220_S25_pair_trim_R2.fastq.gz	$output/W154_pair_trim_R2.fastq.gz
cp	$input2/POS_CTRL_7_S55_pair_trim_R2.fastq.gz	$output/W155_pair_trim_R2.fastq.gz
cp	$input2/COLL_232_S50_pair_trim_R2.fastq.gz	$output/W156_pair_trim_R2.fastq.gz
cp	$input2/COLL_200_S13_pair_trim_R2.fastq.gz	$output/W157_pair_trim_R2.fastq.gz
cp	$input2/COLL_97_S38_pair_trim_R2.fastq.gz	$output/W158_pair_trim_R2.fastq.gz
cp	$input2/COLL_219_S23_pair_trim_R2.fastq.gz	$output/W159_pair_trim_R2.fastq.gz
cp	$input2/COLL_176_S8_pair_trim_R2.fastq.gz	$output/W160_pair_trim_R2.fastq.gz
cp	$input2/COLL_208_S17_pair_trim_R2.fastq.gz	$output/W161_pair_trim_R2.fastq.gz
cp	$input2/COLL_231_S49_pair_trim_R2.fastq.gz	$output/W162_pair_trim_R2.fastq.gz
cp	$input2/COLL_177_S9_pair_trim_R2.fastq.gz	$output/W163_pair_trim_R2.fastq.gz
cp	$input2/COLL_201_S14_pair_trim_R2.fastq.gz	$output/W164_pair_trim_R2.fastq.gz
cp	$input2/COLL_234_S52_pair_trim_R2.fastq.gz	$output/W165_pair_trim_R2.fastq.gz
cp	$input2/NO_147_S33_pair_trim_R2.fastq.gz	$output/W166_pair_trim_R2.fastq.gz
cp	$input2/COLL_205_S16_pair_trim_R2.fastq.gz	$output/W167_pair_trim_R2.fastq.gz
cp	$input2/COLL_179_S10_pair_trim_R2.fastq.gz	$output/W168_pair_trim_R2.fastq.gz
cp	$input2/COLL_174_S6_pair_trim_R2.fastq.gz	$output/W169_pair_trim_R2.fastq.gz
cp	$input2/COLL_160_S1_pair_trim_R2.fastq.gz	$output/W170_pair_trim_R2.fastq.gz
cp	$input2/COLL_173_S5_pair_trim_R2.fastq.gz	$output/W171_pair_trim_R2.fastq.gz
cp	$input2/COLL_226_S44_pair_trim_R2.fastq.gz	$output/W172_pair_trim_R2.fastq.gz
cp	$input2/COLL_57_S41_pair_trim_R2.fastq.gz	$output/W173_pair_trim_R2.fastq.gz
cp	$input2/COLL_218_S22_pair_trim_R2.fastq.gz	$output/W174_pair_trim_R2.fastq.gz
cp	$input2/NEG_CTRL_7_S48_pair_trim_R2.fastq.gz	$output/W175_pair_trim_R2.fastq.gz
cp	$input2/COLL_213_S43_pair_trim_R2.fastq.gz	$output/W176_pair_trim_R2.fastq.gz
cp	$input2/COLL_236_S54_pair_trim_R2.fastq.gz	$output/W177_pair_trim_R2.fastq.gz
cp	$input2/COLL_118_S42_pair_trim_R2.fastq.gz	$output/W178_pair_trim_R2.fastq.gz