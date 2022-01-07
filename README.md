# MappingIntHBRV #
In vivo and in vitro mapping of HBRV integrations on the human genome<br><br>

A series of integration sites of HBRV along the human genome were determined experimentally, in vitro and in vivo, using LM-PCR. To test a hypothesis of randon distribution of integrations in the genome, the same number of integrations were simulated according to a random distribution. The integration frecuency at discrete distances from the closest transcriptional start site (TSS) or from the closest CpG Island was determined and compared against the frequency of experimental integrations. <br>

Results of this project have been published at:

## Description of scripts  <br>

* compare_random_exp_integrations_CpGdist.py: Will map the distance to closest CpG Island of experimental and randomly generated integrations. <br>
* compare_random_exp_integrations_TSSdist.py: Will map the distance to closest TSS of experimental and randomly generated integrations. <br>
* compare_random_exp_integrations_TSSdist_-10K.py: The same as previous one, but plotting is done using max. dist. of -10Kb <br>
* compare_random_exp_integrations_TSSdist__60K.py: The same as previous one, but plotting is done using max. dist. of 60Kb <br>
* sum_chr_known.pl: Scripts required to aggregate data to be plotted with script draw_ideogram.R <br>
* draw_ideogram.R: Script to plot the frequency of integration for a window of size 's'. Input file is generate by previous script<br>

