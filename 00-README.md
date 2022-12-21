# salmon-pollinator-visits
Data and code for the article "Marine-derived nutrients from salmon indirectly affect pollinator visits in a coastal wildflower meadow" by Allison M. Dennert, Elizabeth Elle, and John D. Reynolds
---

This dataset is associated with a long-term field experiment conducted in Haíɫzaqv Territory on the Central Coast of British Columbia, Canada. The experiment has a randomized block design with four treatments, involving adding a marine organism to a 1 x 1 m plot, and 25 replicates per treatment. Floral availability (number of blooming flowers) and pollinator visits (number of landing insects) were observed over monthly (June/July/August) 10-minute surveys across two years (2017-2018). We sought to understand the relationship between marine-derived nutrients from salmon and drift seaweed, floral availability, and pollinator visits.

## Description of the Data and file structure

There are two data files included: plantphenology.csv, and pollinator.csv

Details for: plantphenology.csv
--------------------

* Description: A comma-delimited file containing the results of floral availability surveys each experimental treatment plot. Before each pollinator visit survey, and throughout the summer growing season, the number of inflorescences of every plant species was counted.

* Format(s): .csv

* Missing Data Codes: missing data is coded as "NA"

* Dimensions: 1097 rows x 28 columns

* Variables: 
	* date: the date of the floral availability survey
	* year: The year (YYYY) that the survey was conducted
	* pheno.period: the number assigned to each floral availability observation period (possible values 1-6), some of which overlap with the pollinator surveys
	* pollinator: the number assigned to each pollinator survey (possible values 1-3, "n") conducted concurrently with the floral availability surveys; if there was no overlap, the value reads "n" 
	* block: The experimental block (possible values 1-25) that each treatment plot belongs to
	* plot: The experimental treatment applied to a 1 x 1 m plot; A = pink salmon (Oncorhynchus gorbuscha) carcass, B = drift seaweed (Fucus distichus), C = both salmon carcass and seaweed, D = control
	* columns 7-28:  these column names are four-letter species codes denoting the plant species; columns values indicate the number of blooming flowers each species hads in each plot during the survey (ie, floral availability);  
	acmi = Achillea millefolium,
	anlu = Angelica lucida,
	assu = Aster subspicatus (synonym: Symphyotrichum subspicatum), 
	caly = Carex lynbyei,
	cami = Castilleja miniata, 
	capl = Carex pluriflora,
	clsi = Calytonia siberica,
	copa = Conioselinum pacificum,
	frca = Fritillaria camschatcensis,
	gatr = Galium trifidum,
	glma = Glaux maritima,
	grass = All grass species pooled,
	juba = Juncus balticus,
	madi = Maianthemum dilatatum,
	oesa = Oenanthe sarmentosa,
	plmc = Plantago macrocarpa,
	plmr = Plantago maritima,
	poan = Potentilla anserina,
	sili = Sisyrinchium littorale,
	trar = Trientalis arctica,
	trma = Triglochin maritima,
	trwo = Trifolium workskioldii
		
Details for: pollinator.csv
--------------------

* Description: A comma-delimited file containing the results of pollinator visit surveys each experimental treatment plot. Surveys were conducted of each experimental treatment plot for a period of 10 minutes in each of June, July, and August in 2017 and 2018.

* Format(s): .csv

* Missing Data Codes: missing data is coded as "NA"

* Dimensions: 1498 rows x 22 columns

* Variables: 
	* date: the date of the pollinator visit tsurvey
	* year: The year (YYYY) that the survey was conducted
	* pheno.period: the number assigned to each floral availability observation period (possible values 1-6), some of which overlap with the pollinator surveys
	* pollinator: the number assigned to each pollinator survey (possible values 1-3, "n") conducted concurrently with the floral availability surveys; if there was no overlap, the value reads "n" 
	* block: The experimental block (possible values 1-25) that each treatment plot belongs to
	* plot: The experimental treatment applied to a 1 x 1 m plot; A = pink salmon (Oncorhynchus gorbuscha) carcass, B = drift seaweed (Fucus distichus), C = both salmon carcass and seaweed, D = control
	* time: the time of day (HH:MM) that each block's survey was conducted in 24-hour time
	* temperature: the air temperature in degrees celsius (°C) measured on a Kestrel weather meter at the start of each block's survey
	* wind.kmhr: the wind speed in km/hr measured on a Kestrel weather meter at the start of each block's survey
	* cloud.cover: a visual estimate of the proportion of the sky covered by cloud (possible values 0.0-1.0)
	* observer: the initials of the observer during each pollinator visit survey (possible values AD, LW, NR)
	* plant.species: the plant species with available flowers in each plot during the visit surveys; these values correspond to the species names listed in the metadata for columns 7-28 of the plantphenology.csv file
	* number. stems: the number of available/open flowers of each plant species in each plot during the visit surveys; these data correspond to the values observed in each floral availability survey in the plantphenology.csv file
	* bumblebees: the number of visits by the "bumblebee" morphospecies group to to each listed plant species in the "plant.species" column
	* blowflies:  the number of visits by the "blow fly" morphospecies group to to each listed plant species in the "plant.species" column
	* hoverflies:  the number of visits by the "hover fly" morphospecies group to to each listed plant species in the "plant.species" column
	* other.flies:  the number of visits by the "other fly" morphospecies group to to each listed plant species in the "plant.species" column
	 * wasps:  the number of visits by the "wasp" morphospecies group to to each listed plant species in the "plant.species" column
	 * mason.bees:  the number of visits by the "mason bee" morphospecies group to to each listed plant species in the "plant.species" column
	 * sweat.mining.bees:  the number of visits by the "sweat and mining bee" morphospecies group to to each listed plant species in the "plant.species" column
	 * other:  the number of visits by the "other" morphospecies group (all other insect and non-insect visits) to to each listed plant species in the "plant.species" column
	 * total.visits: the sum of the visits by all of the pollinator morphospecies groups to each listed plant species in the "plant.species" column 

## Sharing/access Information

Links to other publicly accessible locations of the data: https://github.com/adennert/salmon-pollinator-visits

These data were generated by the study's co-authors and not derived from another source.







