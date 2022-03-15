

# Essential genes G4s
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/reallyHighGenes.bed | sort   > ActivelyTranscribedGenes/reallyHighGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/highGenes.bed | sort   > ActivelyTranscribedGenes/highGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/midGenes.bed | sort   > ActivelyTranscribedGenes/midGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/lowGenes.bed | sort   > ActivelyTranscribedGenes/lowGenes_G4_K.bed

bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/plusReallyHighGenes.bed | sort   > ActivelyTranscribedGenes/plusreallyHighGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/plusHighGenes.bed | sort   > ActivelyTranscribedGenes/plushighGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/plusMidGenes.bed | sort   > ActivelyTranscribedGenes/plusmidGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/plusLowGenes.bed | sort   > ActivelyTranscribedGenes/pluslowGenes_G4_K.bed

bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/minusReallyHighGenes.bed | sort   > ActivelyTranscribedGenes/minusreallyHighGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/minusHighGenes.bed | sort   > ActivelyTranscribedGenes/minushighGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/minusMidGenes.bed | sort   > ActivelyTranscribedGenes/minusmidGenes_G4_K.bed
bedtools intersect -wa -u -b G4s/G4_K.bed -a ActivelyTranscribedGenes/minusLowGenes.bed | sort   > ActivelyTranscribedGenes/minuslowGenes_G4_K.bed



bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/reallyHighGenes.bed | sort   > ActivelyTranscribedGenes/reallyHighGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/highGenes.bed | sort   > ActivelyTranscribedGenes/highGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/midGenes.bed | sort   > ActivelyTranscribedGenes/midGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/lowGenes.bed | sort   > ActivelyTranscribedGenes/lowGenes_without_G4_K.bed

bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/plusReallyHighGenes.bed | sort   > ActivelyTranscribedGenes/plusreallyHighGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/plusHighGenes.bed | sort   > ActivelyTranscribedGenes/plushighGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/plusMidGenes.bed | sort   > ActivelyTranscribedGenes/plusmidGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/plusLowGenes.bed | sort   > ActivelyTranscribedGenes/pluslowGenes_without_G4_K.bed

bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/minusReallyHighGenes.bed | sort   > ActivelyTranscribedGenes/minusreallyHighGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/minusHighGenes.bed | sort   > ActivelyTranscribedGenes/minushighGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/minusMidGenes.bed | sort   > ActivelyTranscribedGenes/minusmidGenes_without_G4_K.bed
bedtools intersect -v  -b G4s/G4_K.bed -a ActivelyTranscribedGenes/minusLowGenes.bed | sort   > ActivelyTranscribedGenes/minuslowGenes_without_G4_K.bed
