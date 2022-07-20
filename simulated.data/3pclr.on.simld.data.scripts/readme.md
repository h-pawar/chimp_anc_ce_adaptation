Scripts

- 3pclr_runscript.sh
- neutral_3pclr_runscript_10finalreps.sh # to generate an additional 10 neutral reps
# apply 3P-CLR to neutral simulated data

- 3pclr_smaller_s_selection_runscript.sh
# apply 3P-CLR to positive selection simulated data

To extract the central windows from the simulation reps, run the below:
we extract the central window centred on position 599999 per rep, as the beneficial SNP was introduced in the centre of the genomic region simulated

```
for i in m_polymorphic_subset_nc_3pclroutput_*.txt; do \
	awk -v c=2 -v t=599999 'NR==1{d=$c-t;d=d<0?-d:d;v=$c;next}{m=$c-t;m=m<0?-m:m}m<d{d=m;v=$c}END{print v}' $i > test_$i ; \
	value=$(<test_$i); \
	awk -v site="$value" '$2==site' $i > /SAN/ugi/chimp_internalbranches/test_data/selcoeff_final_1000reps/mutations_out/subset_nc/genpos_morgans/3pclr_output/windows/window_$i; \
	done
```

# eg rep, output of central window - 
head window_m_polymorphic_subset_nc_3pclroutput_885_0.1.txt
23 	600616	0.00576591	471169	730863	0.00452322	0.00701628	337.217	0.01	119.356	0.01	53.0551	0.01	497.621	0.1	242.049	0.01
# cols correspond to: #Chr PhysPos	GenPos	PhysStart PhysEnd	GenStart	GenEnd	3PCLR.Anc	3PCLR.Anc.S	3PCLR.A	3PCLR.A.S	3PCLR.B	3PCLR.B.S	XPCLRAC	XPCLRAC.S	XPCLRBC	XPCLRBC.S
